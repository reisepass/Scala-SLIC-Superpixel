package ch.ethz.dalab.scalaslic

import breeze.linalg._
import breeze.numerics._
import scala.collection.mutable

object ImageSegmentation {
  case class SegmentationParams(
    minSize: Int = 100,
    maxSize: Int = 10000,
    colorWeight: Double = 0.5,
    spatialWeight: Double = 0.5,
    textureWeight: Double = 0.3,
    similarityThreshold: Double = 0.7
  )

  def mergeSuperpixels(
    image: Array[Array[Array[Int]]], 
    superpixelMask: Array[Array[Int]], 
    params: SegmentationParams = SegmentationParams()
  ): Array[Array[Int]] = {
    
    // Get unique superpixel IDs
    val superpixelIds = superpixelMask.flatten.distinct.sorted
    
    // Extract features for each superpixel
    val features = superpixelIds.map { id =>
      id -> SuperpixelFeatures.extractFeatures(image, superpixelMask, id)
    }.toMap
    
    // Build adjacency graph
    val adjacencyGraph = buildAdjacencyGraph(superpixelMask)
    
    // Priority queue for merging
    val mergeQueue = new mutable.PriorityQueue[(Int, Int, Double)]()(Ordering.by(-_._3))
    
    // Initialize merge queue with adjacent superpixel pairs
    for {
      id1 <- superpixelIds
      id2 <- adjacencyGraph(id1)
      if id1 < id2
    } {
      val similarity = SuperpixelFeatures.computeSuperpixelSimilarity(features(id1), features(id2))
      mergeQueue.enqueue((id1, id2, similarity))
    }
    
    // Disjoint set for tracking merged regions
    val disjointSet = new DisjointSet(superpixelIds.max + 1)
    
    // Track sizes of merged regions
    val regionSizes = mutable.Map[Int, Int]()
    superpixelIds.foreach { id =>
      regionSizes(id) = features(id).area
    }
    
    // Merge regions
    while (mergeQueue.nonEmpty) {
      val (id1, id2, similarity) = mergeQueue.dequeue()
      
      val root1 = disjointSet.find(id1)
      val root2 = disjointSet.find(id2)
      
      if (root1 != root2 && similarity > params.similarityThreshold) {
        val mergedSize = regionSizes(root1) + regionSizes(root2)
        
        if (mergedSize <= params.maxSize) {
          disjointSet.union(root1, root2)
          val newRoot = disjointSet.find(root1)
          regionSizes(newRoot) = mergedSize
          
          // Update adjacency relationships
          val neighbors = (adjacencyGraph(root1) ++ adjacencyGraph(root2)).filter(_ != root1).filter(_ != root2)
          for (neighbor <- neighbors) {
            if (disjointSet.find(neighbor) != newRoot) {
              val neighborRoot = disjointSet.find(neighbor)
              val newSimilarity = SuperpixelFeatures.computeSuperpixelSimilarity(
                features(root1), features(neighborRoot)
              )
              mergeQueue.enqueue((newRoot, neighborRoot, newSimilarity))
            }
          }
        }
      }
    }
    
    // Create final segmentation mask
    Array.tabulate(superpixelMask.length, superpixelMask(0).length) { (x, y) =>
      disjointSet.find(superpixelMask(x)(y))
    }
  }
  
  def watershedSegmentation(
    image: Array[Array[Array[Int]]], 
    markers: Array[Array[Int]] = null
  ): Array[Array[Int]] = {
    // Convert to grayscale if needed
    val grayscale = if (image(0)(0).length > 1) {
      Array.tabulate(image.length, image(0).length) { (x, y) =>
        (image(x)(y)(0) * 0.299 + image(x)(y)(1) * 0.587 + image(x)(y)(2) * 0.114).toInt
      }
    } else {
      Array.tabulate(image.length, image(0).length) { (x, y) => image(x)(y)(0) }
    }
    
    // Compute gradient magnitude
    val (gradX, gradY, magnitude) = EdgeDetector.sobelOperator(image)
    
    // Initialize markers if not provided
    val initialMarkers = if (markers == null) {
      val localMinima = findLocalMinima(magnitude)
      Array.tabulate(image.length, image(0).length) { (x, y) =>
        if (localMinima(x)(y)) localMinima.flatten.count(identity) else -1
      }
    } else {
      markers
    }
    
    // Priority queue for flooding
    case class Pixel(x: Int, y: Int, value: Double)
    val queue = new mutable.PriorityQueue[Pixel]()(Ordering.by[Pixel, Double](_.value).reverse)
    
    // Initialize queue with marker boundaries
    val result = Array.tabulate(image.length, image(0).length) { (x, y) => initialMarkers(x)(y) }
    for {
      x <- result.indices
      y <- result(0).indices
      if result(x)(y) >= 0
      dx <- -1 to 1
      dy <- -1 to 1
      nx = x + dx
      ny = y + dy
      if nx >= 0 && nx < result.length && ny >= 0 && ny < result(0).length
      if result(nx)(ny) == -1
    } {
      queue.enqueue(Pixel(nx, ny, magnitude(x)(y)))
    }
    
    // Flooding process
    while (queue.nonEmpty) {
      val Pixel(x, y, _) = queue.dequeue()
      
      if (result(x)(y) == -1) {
        // Find labels of neighboring pixels
        val neighborLabels = for {
          dx <- -1 to 1
          dy <- -1 to 1
          if dx != 0 || dy != 0
          nx = x + dx
          ny = y + dy
          if nx >= 0 && nx < result.length && ny >= 0 && ny < result(0).length
          if result(nx)(ny) >= 0
        } yield result(nx)(ny)
        
        // Assign label if all neighbors have the same label
        if (neighborLabels.nonEmpty && neighborLabels.distinct.size == 1) {
          result(x)(y) = neighborLabels.head
          
          // Add neighbors to queue
          for {
            dx <- -1 to 1
            dy <- -1 to 1
            if dx != 0 || dy != 0
            nx = x + dx
            ny = y + dy
            if nx >= 0 && nx < result.length && ny >= 0 && ny < result(0).length
            if result(nx)(ny) == -1
          } {
            queue.enqueue(Pixel(nx, ny, magnitude(nx)(ny)))
          }
        }
      }
    }
    
    result
  }
  
  private def buildAdjacencyGraph(superpixelMask: Array[Array[Int]]): Map[Int, Set[Int]] = {
    val graph = mutable.Map[Int, mutable.Set[Int]]()
    
    for {
      x <- superpixelMask.indices
      y <- superpixelMask(0).indices
      dx <- -1 to 1
      dy <- -1 to 1
      if dx != 0 || dy != 0
      nx = x + dx
      ny = y + dy
      if nx >= 0 && nx < superpixelMask.length && ny >= 0 && ny < superpixelMask(0).length
      if superpixelMask(x)(y) != superpixelMask(nx)(ny)
    } {
      val id1 = superpixelMask(x)(y)
      val id2 = superpixelMask(nx)(ny)
      graph.getOrElseUpdate(id1, mutable.Set[Int]()).add(id2)
      graph.getOrElseUpdate(id2, mutable.Set[Int]()).add(id1)
    }
    
    graph.map { case (k, v) => k -> v.toSet }.toMap
  }
  
  private def findLocalMinima(image: Array[Array[Double]]): Array[Array[Boolean]] = {
    Array.tabulate(image.length, image(0).length) { (x, y) =>
      val value = image(x)(y)
      val isMinimum = (-1 to 1).forall { dx =>
        (-1 to 1).forall { dy =>
          val nx = x + dx
          val ny = y + dy
          if (nx >= 0 && nx < image.length && ny >= 0 && ny < image(0).length) {
            value <= image(nx)(ny)
          } else true
        }
      }
      isMinimum
    }
  }
  
  private class DisjointSet(size: Int) {
    private val parent = Array.range(0, size)
    private val rank = Array.fill(size)(0)
    
    def find(x: Int): Int = {
      if (parent(x) != x) {
        parent(x) = find(parent(x))
      }
      parent(x)
    }
    
    def union(x: Int, y: Int): Unit = {
      val px = find(x)
      val py = find(y)
      if (px != py) {
        if (rank(px) < rank(py)) {
          parent(px) = py
        } else if (rank(px) > rank(py)) {
          parent(py) = px
        } else {
          parent(py) = px
          rank(px) += 1
        }
      }
    }
  }
} 