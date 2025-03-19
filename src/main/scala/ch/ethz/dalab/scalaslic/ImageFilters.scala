package ch.ethz.dalab.scalaslic

import breeze.numerics._
import scala.collection.parallel.CollectionConverters._

object ImageFilters {
  def gaussianKernel(size: Int, sigma: Double): Array[Array[Double]] = {
    val kernel = Array.ofDim[Double](size, size)
    val center = size / 2
    val twoSigmaSquare = 2.0 * sigma * sigma
    
    for {
      x <- 0 until size
      y <- 0 until size
    } {
      val xDistance = x - center
      val yDistance = y - center
      kernel(x)(y) = exp(-(xDistance * xDistance + yDistance * yDistance) / twoSigmaSquare)
    }
    
    val sum = kernel.flatten.sum
    kernel.map(_.map(_ / sum))
  }

  def gaussianBlur(image: Array[Array[Array[Int]]], kernelSize: Int, sigma: Double): Array[Array[Array[Int]]] = {
    val kernel = gaussianKernel(kernelSize, sigma)
    val padSize = kernelSize / 2
    val paddedImage = padImage(image, padSize)
    
    val result = Array.ofDim[Int](image.length, image(0).length, image(0)(0).length)
    
    for {
      x <- 0 until image.length
      y <- 0 until image(0).length
      z <- 0 until image(0)(0).length
    } {
      var sum = 0.0
      for {
        kx <- 0 until kernelSize
        ky <- 0 until kernelSize
      } {
        sum += paddedImage(x + kx)(y + ky)(z) * kernel(kx)(ky)
      }
      result(x)(y)(z) = sum.toInt
    }
    result
  }

  def medianFilter(image: Array[Array[Array[Int]]], windowSize: Int): Array[Array[Array[Int]]] = {
    val padSize = windowSize / 2
    val paddedImage = padImage(image, padSize)
    val result = Array.ofDim[Int](image.length, image(0).length, image(0)(0).length)
    
    for {
      x <- 0 until image.length
      y <- 0 until image(0).length
      z <- 0 until image(0)(0).length
    } {
      val window = for {
        wx <- 0 until windowSize
        wy <- 0 until windowSize
      } yield paddedImage(x + wx)(y + wy)(z)
      
      result(x)(y)(z) = window.sorted.apply(window.size / 2)
    }
    result
  }

  def bilateralFilter(image: Array[Array[Array[Int]]], 
                     windowSize: Int, 
                     sigmaSpace: Double, 
                     sigmaColor: Double): Array[Array[Array[Int]]] = {
    val padSize = windowSize / 2
    val paddedImage = padImage(image, padSize)
    val result = Array.ofDim[Int](image.length, image(0).length, image(0)(0).length)
    
    val gaussian = gaussianKernel(windowSize, sigmaSpace)
    val twoSigmaColorSquare = 2.0 * sigmaColor * sigmaColor
    
    for {
      x <- 0 until image.length
      y <- 0 until image(0).length
      z <- 0 until image(0)(0).length
    } {
      var sum = 0.0
      var weightSum = 0.0
      
      for {
        wx <- 0 until windowSize
        wy <- 0 until windowSize
      } {
        val spatialWeight = gaussian(wx)(wy)
        val colorDiff = paddedImage(x + wx)(y + wy)(z) - paddedImage(x + padSize)(y + padSize)(z)
        val colorWeight = exp(-(colorDiff * colorDiff) / twoSigmaColorSquare)
        val weight = spatialWeight * colorWeight
        
        sum += paddedImage(x + wx)(y + wy)(z) * weight
        weightSum += weight
      }
      
      result(x)(y)(z) = (sum / weightSum).toInt
    }
    result
  }

  private def padImage(image: Array[Array[Array[Int]]], padSize: Int): Array[Array[Array[Int]]] = {
    val paddedImage = Array.ofDim[Int](
      image.length + 2 * padSize,
      image(0).length + 2 * padSize,
      image(0)(0).length
    )
    
    // Copy original image
    for {
      x <- 0 until image.length
      y <- 0 until image(0).length
      z <- 0 until image(0)(0).length
    } {
      paddedImage(x + padSize)(y + padSize)(z) = image(x)(y)(z)
    }
    
    // Pad borders by replicating edge values
    for {
      x <- 0 until paddedImage.length
      y <- 0 until paddedImage(0).length
      z <- 0 until paddedImage(0)(0).length
    } {
      val sourceX = math.min(math.max(x - padSize, 0), image.length - 1)
      val sourceY = math.min(math.max(y - padSize, 0), image(0).length - 1)
      paddedImage(x)(y)(z) = image(sourceX)(sourceY)(z)
    }
    
    paddedImage
  }
} 