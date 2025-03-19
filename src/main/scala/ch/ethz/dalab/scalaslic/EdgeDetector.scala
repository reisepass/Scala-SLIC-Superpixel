package ch.ethz.dalab.scalaslic

import breeze.numerics._

object EdgeDetector {
  val sobelX = Array(
    Array(-1.0, 0.0, 1.0),
    Array(-2.0, 0.0, 2.0),
    Array(-1.0, 0.0, 1.0)
  )

  val sobelY = Array(
    Array(-1.0, -2.0, -1.0),
    Array(0.0, 0.0, 0.0),
    Array(1.0, 2.0, 1.0)
  )

  def sobelOperator(image: Array[Array[Array[Int]]]): (Array[Array[Double]], Array[Array[Double]], Array[Array[Double]]) = {
    val width = image.length
    val height = image(0).length
    val gradientX = Array.ofDim[Double](width, height)
    val gradientY = Array.ofDim[Double](width, height)
    val magnitude = Array.ofDim[Double](width, height)

    for {
      x <- 1 until width - 1
      y <- 1 until height - 1
    } {
      var sumX = 0.0
      var sumY = 0.0
      
      for {
        i <- -1 to 1
        j <- -1 to 1
      } {
        val pixel = image(x + i)(y + j)(0) // Using grayscale or first channel
        sumX += pixel * sobelX(i + 1)(j + 1)
        sumY += pixel * sobelY(i + 1)(j + 1)
      }
      
      gradientX(x)(y) = sumX
      gradientY(x)(y) = sumY
      magnitude(x)(y) = sqrt(sumX * sumX + sumY * sumY)
    }

    (gradientX, gradientY, magnitude)
  }

  def cannyEdgeDetector(image: Array[Array[Array[Int]]], 
                       lowThreshold: Double = 20, 
                       highThreshold: Double = 50): Array[Array[Boolean]] = {
    // Step 1: Apply Gaussian blur to reduce noise
    val blurred = ImageFilters.gaussianBlur(image, 5, 1.4)
    
    // Step 2: Compute gradients using Sobel
    val (gradX, gradY, magnitude) = sobelOperator(blurred)
    
    // Step 3: Non-maximum suppression
    val width = image.length
    val height = image(0).length
    val suppressed = Array.ofDim[Double](width, height)
    val edges = Array.ofDim[Boolean](width, height)
    
    for {
      x <- 1 until width - 1
      y <- 1 until height - 1
    } {
      val angle = atan2(gradY(x)(y), gradX(x)(y)) * 180 / math.Pi
      val dir = ((angle + 180) / 45).toInt % 4
      
      val mag = magnitude(x)(y)
      var isMax = false
      
      // Check if pixel is local maximum in gradient direction
      dir match {
        case 0 => isMax = mag > magnitude(x)(y-1) && mag > magnitude(x)(y+1)
        case 1 => isMax = mag > magnitude(x-1)(y-1) && mag > magnitude(x+1)(y+1)
        case 2 => isMax = mag > magnitude(x-1)(y) && mag > magnitude(x+1)(y)
        case 3 => isMax = mag > magnitude(x-1)(y+1) && mag > magnitude(x+1)(y-1)
      }
      
      suppressed(x)(y) = if (isMax) mag else 0
    }
    
    // Step 4: Double threshold and hysteresis
    for {
      x <- 1 until width - 1
      y <- 1 until height - 1
    } {
      if (suppressed(x)(y) >= highThreshold) {
        edges(x)(y) = true
        // Trace edges through pixels above low threshold
        def trace(i: Int, j: Int): Unit = {
          for {
            di <- -1 to 1
            dj <- -1 to 1
            if di != 0 || dj != 0
          } {
            val ni = i + di
            val nj = j + dj
            if (ni >= 0 && ni < width && nj >= 0 && nj < height &&
                !edges(ni)(nj) && suppressed(ni)(nj) >= lowThreshold) {
              edges(ni)(nj) = true
              trace(ni, nj)
            }
          }
        }
        trace(x, y)
      }
    }
    
    edges
  }
} 