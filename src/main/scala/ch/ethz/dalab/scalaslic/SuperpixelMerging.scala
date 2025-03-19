package ch.ethz.dalab.scalaslic

import scala.collection.mutable
import breeze.linalg._
import breeze.numerics._

object SuperpixelMerging {
  case class MergeParams(
    minSize: Int = 50,
    maxSize: Int = 1000,
    colorSimilarityWeight: Double = 0.6,
    textureSimilarityWeight: Double = 0.2,
    spatialSimilarityWeight: Double = 0.2,
    similarityThreshold: Double = 0.7
  )

  def mergeSmallSuperpixels(
    image: Array[Array[Array[Int]]], 
    superpixelMask: Array[Array[Int]], 
    params: MergeParams = MergeParams()
  ): Array[Array[Int]] = {
    
    // Get superpixel sizes
    val sizes = mutable.Map[Int, Int]()
    for {
      x <- superpixelMask.indices
      y <- superpixelMask(0).indices
    } {
      val id = superpixelMask(x)(y)
      sizes(id) = sizes.getOrElse(id, 0) + 1
    }
    
    // Find small superpixels
    val smallSuperpixels = sizes.filter(_._2 < params.minSize).keys.toSet
    
    if (smallSuperpixels.isEmpty) {
      return superpixelMask
    }
    
    // Build adjacency graph
    val adjacencyGraph = buildAdjacencyGraph(superpixelMask)
    
    // Extract features for all superpixels
    val features = superpixelMask.flatten.distinct.map { id =>
      id -> SuperpixelFeatures.extractFeatures(image, superpixelMask, id)
    }.toMap
    
    // Process small superpixels
    val result = Array.tabulate(superpixelMask.length, superpixelMask(0).length) { (x, y) =>
      superpixelMask(x)(y)
    }
    
    for (smallId <- smallSuperpixels) {
      if (sizes(smallId) < params.minSize) {
        // Find best neighbor to merge with
        val neighbors = adjacencyGraph.getOrElse(smallId, Set.empty)
        val validNeighbors = neighbors.filter(n => !smallSuperpixels.contains(n) || sizes(n) + sizes(smallId) <= params.maxSize)
        
        if (validNeighbors.nonEmpty) {
          val bestNeighbor = validNeighbors.maxBy { neighborId =>
            computeMergeSimilarity(
              features(smallId), 
              features(neighborId),
              params
            )
          }
          
          // Perform merge
          for {
            x <- result.indices
            y <- result(0).indices
            if result(x)(y) == smallId
          } {
            result(x)(y) = bestNeighbor
          }
          
          // Update sizes
          sizes(bestNeighbor) += sizes(smallId)
          sizes.remove(smallId)
        }
      }
    }
    
    result
  }

  def enforceConnectivity(
    superpixelMask: Array[Array[Int]], 
    minSize: Int = 20
  ): Array[Array[Int]] = {
    val width = superpixelMask.length
    val height = superpixelMask(0).length
    val result = Array.tabulate(width, height) { (x, y) => superpixelMask(x)(y) }
    val visited = Array.ofDim[Boolean](width, height)
    var nextLabel = superpixelMask.flatten.max + 1
    
    def floodFill(startX: Int, startY: Int, label: Int): Int = {
      val queue = mutable.Queue[(Int, Int)]()
      queue.enqueue((startX, startY))
      var size = 0
      
      while (queue.nonEmpty) {
        val (x, y) = queue.dequeue()
        if (!visited(x)(y) && result(x)(y) == label) {
          visited(x)(y) = true
          size += 1
          
          // Add neighbors
          for {
            dx <- -1 to 1
            dy <- -1 to 1
            if dx != 0 || dy != 0
            nx = x + dx
            ny = y + dy
            if nx >= 0 && nx < width && ny >= 0 && ny < height
          } {
            if (!visited(nx)(ny) && result(nx)(ny) == label) {
              queue.enqueue((nx, ny))
            }
          }
        }
      }
      size
    }
    
    // Process each pixel
    for {
      x <- result.indices
      y <- result(0).indices
      if !visited(x)(y)
    } {
      val label = result(x)(y)
      val size = floodFill(x, y, label)
      
      // If component is too small, merge with largest neighbor
      if (size < minSize) {
        val neighborLabels = mutable.Map[Int, Int]()
        
        // Count neighbor labels
        for {
          dx <- -1 to 1
          dy <- -1 to 1
          if dx != 0 || dy != 0
          nx = x + dx
          ny = y + dy
          if nx >= 0 && nx < width && ny >= 0 && ny < height
          neighborLabel = result(nx)(ny)
          if neighborLabel != label
        } {
          neighborLabels(neighborLabel) = neighborLabels.getOrElse(neighborLabel, 0) + 1
        }
        
        if (neighborLabels.nonEmpty) {
          val newLabel = neighborLabels.maxBy(_._2)._1
          
          // Relabel small component
          for {
            px <- result.indices
            py <- result(0).indices
            if result(px)(py) == label
          } {
            result(px)(py) = newLabel
          }
        } else {
          // Isolated component, assign new label
          for {
            px <- result.indices
            py <- result(0).indices
            if result(px)(py) == label
          } {
            result(px)(py) = nextLabel
          }
          nextLabel += 1
        }
      }
    }
    
    result
  }

  def smoothBoundaries(
    image: Array[Array[Array[Int]]], 
    superpixelMask: Array[Array[Int]], 
    windowSize: Int = 3
  ): Array[Array[Int]] = {
    val result = Array.tabulate(superpixelMask.length, superpixelMask(0).length) { (x, y) =>
      superpixelMask(x)(y)
    }
    val radius = windowSize / 2
    
    // Find boundary pixels
    val boundaryPixels = for {
      x <- radius until superpixelMask.length - radius
      y <- radius until superpixelMask(0).length - radius
      if isBoundaryPixel(superpixelMask, x, y)
    } yield (x, y)
    
    // Process boundary pixels
    for ((x, y) <- boundaryPixels) {
      val currentLabel = result(x)(y)
      val neighborLabels = mutable.Map[Int, Double]()
      
      // Calculate color-weighted votes for each label in the window
      for {
        wx <- -radius to radius
        wy <- -radius to radius
        nx = x + wx
        ny = y + wy
      } {
        val neighborLabel = result(nx)(ny)
        if (neighborLabel != currentLabel) {
          val colorDiff = (0 until image(0)(0).length).map { c =>
            math.pow(image(x)(y)(c) - image(nx)(ny)(c), 2)
          }.sum
          val weight = math.exp(-colorDiff / 100.0)
          neighborLabels(neighborLabel) = neighborLabels.getOrElse(neighborLabel, 0.0) + weight
        }
      }
      
      // Assign most similar label
      if (neighborLabels.nonEmpty) {
        val newLabel = neighborLabels.maxBy(_._2)._1
        result(x)(y) = newLabel
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

  private def computeMergeSimilarity(
    features1: SuperpixelFeatures,
    features2: SuperpixelFeatures,
    params: MergeParams
  ): Double = {
    // Color similarity
    val colorDiff = (features1.meanColor, features2.meanColor).zipped.map((a, b) => (a - b) * (a - b)).sum
    val colorSim = exp(-colorDiff / 100.0)
    
    // Texture similarity using gradient histograms
    val textureDiff = (features1.textureGradients, features2.textureGradients).zipped.map((a, b) => (a - b) * (a - b)).sum
    val textureSim = exp(-textureDiff / 1000.0)
    
    // Spatial proximity
    val spatialDiff = math.pow(features1.centroid._1 - features2.centroid._1, 2) + 
                     math.pow(features1.centroid._2 - features2.centroid._2, 2)
    val spatialSim = exp(-spatialDiff / (100.0 * 100.0))
    
    // Weighted combination
    params.colorSimilarityWeight * colorSim +
    params.textureSimilarityWeight * textureSim +
    params.spatialSimilarityWeight * spatialSim
  }

  private def isBoundaryPixel(mask: Array[Array[Int]], x: Int, y: Int): Boolean = {
    val label = mask(x)(y)
    (-1 to 1).exists { dx =>
      (-1 to 1).exists { dy =>
        if (dx != 0 || dy != 0) {
          val nx = x + dx
          val ny = y + dy
          nx >= 0 && nx < mask.length && ny >= 0 && ny < mask(0).length && mask(nx)(ny) != label
        } else false
      }
    }
  }
} 