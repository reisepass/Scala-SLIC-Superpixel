package ch.ethz.dalab.scalaslic

import breeze.numerics._
import scala.collection.parallel.CollectionConverters._

object ImageEnhancement {
  def histogramEqualization(image: Array[Array[Array[Int]]]): Array[Array[Array[Int]]] = {
    val result = Array.ofDim[Int](image.length, image(0).length, image(0)(0).length)
    
    // Process each channel independently
    for (channel <- image(0)(0).indices) {
      // Calculate histogram
      val histogram = Array.fill(256)(0)
      for {
        x <- image.indices
        y <- image(0).indices
      } {
        histogram(image(x)(y)(channel)) += 1
      }
      
      // Calculate cumulative distribution function
      val cdf = histogram.scanLeft(0)(_ + _).tail
      val cdfMin = cdf.filter(_ > 0).min
      val denominator = image.length * image(0).length - cdfMin
      
      // Apply equalization
      for {
        x <- image.indices
        y <- image(0).indices
      } {
        result(x)(y)(channel) = ((cdf(image(x)(y)(channel)) - cdfMin) * 255.0 / denominator).round.toInt
      }
    }
    
    result
  }

  def adaptiveHistogramEqualization(
    image: Array[Array[Array[Int]]], 
    tileSize: Int = 8, 
    clipLimit: Double = 4.0
  ): Array[Array[Array[Int]]] = {
    val result = Array.ofDim[Int](image.length, image(0).length, image(0)(0).length)
    
    // Process each channel independently
    for (channel <- image(0)(0).indices) {
      val tiles = for {
        tx <- 0 until (image.length + tileSize - 1) / tileSize
        ty <- 0 until (image(0).length + tileSize - 1) / tileSize
      } yield {
        val startX = tx * tileSize
        val startY = ty * tileSize
        val endX = math.min(startX + tileSize, image.length)
        val endY = math.min(startY + tileSize, image(0).length)
        
        // Calculate histogram for tile
        val histogram = Array.fill(256)(0)
        for {
          x <- startX until endX
          y <- startY until endY
        } {
          histogram(image(x)(y)(channel)) += 1
        }
        
        // Clip histogram
        val clipHeight = (clipLimit * (endX - startX) * (endY - startY) / 256.0).toInt
        val excess = histogram.map(h => math.max(0, h - clipHeight)).sum
        val redistribution = excess / 256
        for (i <- histogram.indices) {
          histogram(i) = math.min(histogram(i), clipHeight) + redistribution
        }
        
        // Calculate mapping function
        val cdf = histogram.scanLeft(0)(_ + _).tail
        val scale = 255.0 / cdf.last
        val mapping = cdf.map(c => (c * scale).toInt)
        
        (startX, startY, endX, endY, mapping)
      }
      
      // Apply mappings with bilinear interpolation
      for {
        x <- image.indices
        y <- image(0).indices
      } {
        val tileX = x / tileSize
        val tileY = y / tileSize
        
        // Find surrounding tiles
        val surroundingTiles = for {
          dx <- 0 to 1
          dy <- 0 to 1
          tx = math.min(tileX + dx, (image.length - 1) / tileSize)
          ty = math.min(tileY + dy, (image(0).length - 1) / tileSize)
          tile <- tiles.find(t => t._1 / tileSize == tx && t._2 / tileSize == ty)
        } yield tile
        
        if (surroundingTiles.length == 1) {
          // Use single tile mapping
          val tile = surroundingTiles.head
          result(x)(y)(channel) = tile._5(image(x)(y)(channel))
        } else {
          // Interpolate between tiles
          val xWeight = (x % tileSize).toDouble / tileSize
          val yWeight = (y % tileSize).toDouble / tileSize
          
          val interpolated = surroundingTiles.map { tile =>
            val tileXWeight = if (tile._1 / tileSize == tileX) 1 - xWeight else xWeight
            val tileYWeight = if (tile._2 / tileSize == tileY) 1 - yWeight else yWeight
            tile._5(image(x)(y)(channel)) * tileXWeight * tileYWeight
          }.sum
          
          result(x)(y)(channel) = interpolated.toInt
        }
      }
    }
    
    result
  }

  def contrastStretching(
    image: Array[Array[Array[Int]]], 
    lowPercentile: Double = 1.0, 
    highPercentile: Double = 99.0
  ): Array[Array[Array[Int]]] = {
    val result = Array.ofDim[Int](image.length, image(0).length, image(0)(0).length)
    
    for (channel <- image(0)(0).indices) {
      // Get channel values and sort them
      val values = (for {
        x <- image.indices
        y <- image(0).indices
      } yield image(x)(y)(channel)).sorted
      
      // Calculate percentile values
      val lowIndex = ((values.length * lowPercentile) / 100.0).toInt
      val highIndex = ((values.length * highPercentile) / 100.0).toInt
      val low = values(lowIndex)
      val high = values(highIndex)
      val range = high - low
      
      // Apply contrast stretching
      for {
        x <- image.indices
        y <- image(0).indices
      } {
        val stretched = if (range > 0) {
          ((image(x)(y)(channel) - low) * 255.0 / range).round.toInt
        } else {
          image(x)(y)(channel)
        }
        result(x)(y)(channel) = math.min(255, math.max(0, stretched))
      }
    }
    
    result
  }

  def unsharpMasking(
    image: Array[Array[Array[Int]]], 
    sigma: Double = 1.0, 
    amount: Double = 1.5
  ): Array[Array[Array[Int]]] = {
    // Create blurred version
    val blurred = ImageFilters.gaussianBlur(image, 5, sigma)
    
    // Apply unsharp mask
    Array.tabulate(image.length, image(0).length, image(0)(0).length) { (x, y, c) =>
      val detail = image(x)(y)(c) - blurred(x)(y)(c)
      val enhanced = image(x)(y)(c) + (detail * amount).toInt
      math.min(255, math.max(0, enhanced))
    }
  }

  def gammaCorrection(image: Array[Array[Array[Int]]], gamma: Double): Array[Array[Array[Int]]] = {
    val lookupTable = Array.tabulate(256)(i => (255 * math.pow(i / 255.0, gamma)).toInt)
    
    Array.tabulate(image.length, image(0).length, image(0)(0).length) { (x, y, c) =>
      lookupTable(image(x)(y)(c))
    }
  }

  def localContrastEnhancement(
    image: Array[Array[Array[Int]]], 
    windowSize: Int = 7, 
    k: Double = 0.4
  ): Array[Array[Array[Int]]] = {
    val result = Array.ofDim[Int](image.length, image(0).length, image(0)(0).length)
    val radius = windowSize / 2
    
    for {
      channel <- image(0)(0).indices
      x <- image.indices
      y <- image(0).indices
    } {
      // Calculate local statistics
      var sum = 0.0
      var sumSq = 0.0
      var count = 0
      
      for {
        wx <- math.max(0, x - radius) until math.min(image.length, x + radius + 1)
        wy <- math.max(0, y - radius) until math.min(image(0).length, y + radius + 1)
      } {
        val value = image(wx)(wy)(channel)
        sum += value
        sumSq += value * value
        count += 1
      }
      
      val mean = sum / count
      val variance = (sumSq / count) - (mean * mean)
      val stdDev = math.sqrt(variance)
      
      // Enhance contrast
      val pixelValue = image(x)(y)(channel)
      val enhanced = mean + k * (pixelValue - mean) * (1 + stdDev)
      result(x)(y)(channel) = math.min(255, math.max(0, enhanced.toInt))
    }
    
    result
  }
} 