package ch.ethz.dalab.scalaslic

import scala.collection.mutable
import breeze.numerics._

object RegionGrowing {
  case class RegionGrowingParams(
    colorThreshold: Double = 30.0,
    minRegionSize: Int = 100,
    maxRegionSize: Int = 1000,
    connectivityType: Int = 8  // 4 or 8 connectivity
  )

  def growRegions(
    image: Array[Array[Array[Int]]], 
    seeds: Array[(Int, Int)] = null,
    params: RegionGrowingParams = RegionGrowingParams()
  ): Array[Array[Int]] = {
    val width = image.length
    val height = image(0).length
    val result = Array.fill(width, height)(-1)
    val visited = Array.fill(width, height)(false)
    
    // Generate seeds if not provided
    val initialSeeds = if (seeds == null) {
      generateSeeds(image, width * height / params.maxRegionSize)
    } else seeds
    
    // Process each seed
    var currentLabel = 0
    for ((seedX, seedY) <- initialSeeds if !visited(seedX)(seedY)) {
      val region = growRegion(
        image, result, visited, seedX, seedY, currentLabel, params
      )
      
      if (region.size >= params.minRegionSize) {
        // Mark region with current label
        for ((x, y) <- region) {
          result(x)(y) = currentLabel
        }
        currentLabel += 1
      }
    }
    
    // Handle unassigned pixels
    handleUnassignedPixels(image, result, params)
    
    result
  }

  def splitAndMerge(
    image: Array[Array[Array[Int]]], 
    minBlockSize: Int = 16,
    similarityThreshold: Double = 30.0
  ): Array[Array[Int]] = {
    val width = image.length
    val height = image(0).length
    val result = Array.fill(width, height)(-1)
    var currentLabel = 0
    
    // Initial splitting phase
    val blocks = splitPhase(image, minBlockSize, similarityThreshold)
    
    // Label initial blocks
    for (block <- blocks) {
      for {
        x <- block.x until block.x + block.width
        y <- block.y until block.y + block.height
      } {
        result(x)(y) = currentLabel
      }
      currentLabel += 1
    }
    
    // Merging phase
    val finalResult = mergePhase(image, result, blocks, similarityThreshold)
    
    // Relabel regions to ensure consecutive labels
    relabelRegions(finalResult)
  }

  private case class Block(
    x: Int, y: Int, 
    width: Int, height: Int, 
    meanColor: Array[Double],
    variance: Double
  )

  private def splitPhase(
    image: Array[Array[Array[Int]]], 
    minBlockSize: Int,
    similarityThreshold: Double
  ): List[Block] = {
    def split(x: Int, y: Int, width: Int, height: Int): List[Block] = {
      if (width <= minBlockSize || height <= minBlockSize) {
        // Calculate block statistics
        val stats = calculateBlockStatistics(image, x, y, width, height)
        List(Block(x, y, width, height, stats._1, stats._2))
      } else {
        val variance = calculateBlockVariance(image, x, y, width, height)
        
        if (variance < similarityThreshold) {
          val stats = calculateBlockStatistics(image, x, y, width, height)
          List(Block(x, y, width, height, stats._1, stats._2))
        } else {
          // Split block into four
          val halfWidth = width / 2
          val halfHeight = height / 2
          
          split(x, y, halfWidth, halfHeight) :::
          split(x + halfWidth, y, width - halfWidth, halfHeight) :::
          split(x, y + halfHeight, halfWidth, height - halfHeight) :::
          split(x + halfWidth, y + halfHeight, width - halfWidth, height - halfHeight)
        }
      }
    }
    
    split(0, 0, image.length, image(0).length)
  }

  private def mergePhase(
    image: Array[Array[Array[Int]]], 
    labels: Array[Array[Int]],
    blocks: List[Block],
    similarityThreshold: Double
  ): Array[Array[Int]] = {
    val result = Array.tabulate(labels.length, labels(0).length) { (x, y) => labels(x)(y) }
    val merged = mutable.Set[(Int, Int)]()
    var changed = true
    
    while (changed) {
      changed = false
      
      for {
        i <- blocks.indices
        j <- blocks.indices
        if i < j && !merged.contains((i, j))
      } {
        val block1 = blocks(i)
        val block2 = blocks(j)
        
        // Check if blocks are adjacent
        if (areBlocksAdjacent(block1, block2)) {
          val similarity = calculateBlockSimilarity(block1, block2)
          
          if (similarity < similarityThreshold) {
            // Merge blocks
            for {
              x <- block2.x until block2.x + block2.width
              y <- block2.y until block2.y + block2.height
              if x < result.length && y < result(0).length
            } {
              result(x)(y) = result(block1.x)(block1.y)
            }
            
            merged.add((i, j))
            changed = true
          }
        }
      }
    }
    
    result
  }

  private def generateSeeds(
    image: Array[Array[Array[Int]]], 
    numSeeds: Int
  ): Array[(Int, Int)] = {
    val width = image.length
    val height = image(0).length
    val gradientMagnitude = calculateGradientMagnitude(image)
    
    // Use gradient information to avoid placing seeds on edges
    val candidates = for {
      x <- 1 until width - 1
      y <- 1 until height - 1
      if gradientMagnitude(x)(y) < 50  // Threshold for edge detection
    } yield (x, y)
    
    // Randomly select seeds from candidates
    val random = new scala.util.Random()
    random.shuffle(candidates.toList).take(numSeeds).toArray
  }

  private def growRegion(
    image: Array[Array[Array[Int]]], 
    labels: Array[Array[Int]],
    visited: Array[Array[Boolean]],
    seedX: Int,
    seedY: Int,
    label: Int,
    params: RegionGrowingParams
  ): Set[(Int, Int)] = {
    val region = mutable.Set[(Int, Int)]()
    val queue = mutable.Queue[(Int, Int)]()
    val seedColor = image(seedX)(seedY)
    
    queue.enqueue((seedX, seedY))
    visited(seedX)(seedY) = true
    region.add((seedX, seedY))
    
    while (queue.nonEmpty && region.size < params.maxRegionSize) {
      val (x, y) = queue.dequeue()
      
      // Get neighbors based on connectivity type
      val neighbors = params.connectivityType match {
        case 4 => List((-1,0), (1,0), (0,-1), (0,1))
        case 8 => List((-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1))
        case _ => throw new IllegalArgumentException("Invalid connectivity type")
      }
      
      for {
        (dx, dy) <- neighbors
        nx = x + dx
        ny = y + dy
        if nx >= 0 && nx < image.length && 
           ny >= 0 && ny < image(0).length && 
           !visited(nx)(ny)
      } {
        val colorDiff = (0 until image(0)(0).length).map { c =>
          math.pow(image(nx)(ny)(c) - seedColor(c), 2)
        }.sum
        
        if (colorDiff < params.colorThreshold) {
          queue.enqueue((nx, ny))
          visited(nx)(ny) = true
          region.add((nx, ny))
        }
      }
    }
    
    region.toSet
  }

  private def handleUnassignedPixels(
    image: Array[Array[Array[Int]]], 
    labels: Array[Array[Int]],
    params: RegionGrowingParams
  ): Unit = {
    val width = image.length
    val height = image(0).length
    
    for {
      x <- 0 until width
      y <- 0 until height
      if labels(x)(y) == -1
    } {
      // Find best matching neighbor region
      var bestLabel = -1
      var minDiff = Double.MaxValue
      
      val neighbors = params.connectivityType match {
        case 4 => List((-1,0), (1,0), (0,-1), (0,1))
        case 8 => List((-1,-1), (-1,0), (-1,1), (0,-1), (0,1), (1,-1), (1,0), (1,1))
        case _ => throw new IllegalArgumentException("Invalid connectivity type")
      }
      
      for {
        (dx, dy) <- neighbors
        nx = x + dx
        ny = y + dy
        if nx >= 0 && nx < width && ny >= 0 && ny < height && labels(nx)(ny) != -1
      } {
        val colorDiff = (0 until image(0)(0).length).map { c =>
          math.pow(image(x)(y)(c) - image(nx)(ny)(c), 2)
        }.sum
        
        if (colorDiff < minDiff) {
          minDiff = colorDiff
          bestLabel = labels(nx)(ny)
        }
      }
      
      if (bestLabel != -1) {
        labels(x)(y) = bestLabel
      }
    }
  }

  private def calculateBlockStatistics(
    image: Array[Array[Array[Int]]], 
    x: Int, y: Int, 
    width: Int, height: Int
  ): (Array[Double], Double) = {
    val meanColor = Array.fill(image(0)(0).length)(0.0)
    var sumSquaredDiff = 0.0
    
    // Calculate mean color
    for {
      i <- x until x + width
      j <- y until y + height
      c <- 0 until image(0)(0).length
    } {
      meanColor(c) += image(i)(j)(c)
    }
    
    val numPixels = width * height
    for (c <- meanColor.indices) {
      meanColor(c) /= numPixels
    }
    
    // Calculate variance
    for {
      i <- x until x + width
      j <- y until y + height
      c <- 0 until image(0)(0).length
    } {
      val diff = image(i)(j)(c) - meanColor(c)
      sumSquaredDiff += diff * diff
    }
    
    val variance = sumSquaredDiff / (numPixels * image(0)(0).length)
    (meanColor, variance)
  }

  private def calculateBlockVariance(
    image: Array[Array[Array[Int]]], 
    x: Int, y: Int, 
    width: Int, height: Int
  ): Double = {
    calculateBlockStatistics(image, x, y, width, height)._2
  }

  private def areBlocksAdjacent(block1: Block, block2: Block): Boolean = {
    val x1Range = block1.x until (block1.x + block1.width)
    val y1Range = block1.y until (block1.y + block1.height)
    val x2Range = block2.x until (block2.x + block2.width)
    val y2Range = block2.y until (block2.y + block2.height)
    
    (x1Range.last + 1 == x2Range.head || x2Range.last + 1 == x1Range.head) && 
    y1Range.intersect(y2Range).nonEmpty ||
    (y1Range.last + 1 == y2Range.head || y2Range.last + 1 == y1Range.head) && 
    x1Range.intersect(x2Range).nonEmpty
  }

  private def calculateBlockSimilarity(block1: Block, block2: Block): Double = {
    val colorDiff = (block1.meanColor, block2.meanColor).zipped.map { (c1, c2) =>
      math.pow(c1 - c2, 2)
    }.sum
    
    val varianceDiff = math.abs(block1.variance - block2.variance)
    colorDiff + varianceDiff
  }

  private def calculateGradientMagnitude(image: Array[Array[Array[Int]]]): Array[Array[Double]] = {
    val width = image.length
    val height = image(0).length
    val result = Array.ofDim[Double](width, height)
    
    for {
      x <- 1 until width - 1
      y <- 1 until height - 1
    } {
      var sumSquaredGradient = 0.0
      
      for (c <- 0 until image(0)(0).length) {
        val gx = image(x+1)(y)(c) - image(x-1)(y)(c)
        val gy = image(x)(y+1)(c) - image(x)(y-1)(c)
        sumSquaredGradient += gx * gx + gy * gy
      }
      
      result(x)(y) = math.sqrt(sumSquaredGradient)
    }
    
    result
  }

  private def relabelRegions(labels: Array[Array[Int]]): Array[Array[Int]] = {
    val oldToNew = labels.flatten.distinct.sorted.zipWithIndex.toMap
    Array.tabulate(labels.length, labels(0).length) { (x, y) =>
      oldToNew(labels(x)(y))
    }
  }
} 