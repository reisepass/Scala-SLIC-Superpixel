package ch.ethz.dalab.scalaslic

import breeze.linalg._
import breeze.numerics._
import breeze.stats._

object TextureAnalysis {
  case class GLCMFeatures(
    contrast: Double,
    correlation: Double,
    energy: Double,
    homogeneity: Double,
    entropy: Double
  )

  def computeGLCM(
    image: Array[Array[Int]], 
    distance: Int = 1, 
    angle: Double = 0.0, 
    levels: Int = 256
  ): DenseMatrix[Double] = {
    val glcm = DenseMatrix.zeros[Double](levels, levels)
    val dx = (distance * math.cos(angle)).round.toInt
    val dy = (distance * math.sin(angle)).round.toInt
    
    for {
      x <- 0 until image.length
      y <- 0 until image(0).length
      nx = x + dx
      ny = y + dy
      if nx >= 0 && nx < image.length && ny >= 0 && ny < image(0).length
    } {
      val i = image(x)(y)
      val j = image(nx)(ny)
      glcm(i, j) += 1
    }
    
    // Normalize
    glcm :/= sum(glcm)
    glcm
  }

  def extractGLCMFeatures(glcm: DenseMatrix[Double]): GLCMFeatures = {
    val size = glcm.rows
    var contrast = 0.0
    var correlation = 0.0
    var energy = 0.0
    var homogeneity = 0.0
    var entropy = 0.0
    
    // Calculate means and standard deviations
    val meanI = sum(DenseVector.tabulate(size)(i => 
      sum(glcm(i, ::).t :* DenseVector.ones[Double](size)) * i
    ))
    val meanJ = sum(DenseVector.tabulate(size)(j => 
      sum(glcm(::, j) :* DenseVector.ones[Double](size)) * j
    ))
    
    val stdI = sqrt(sum(DenseVector.tabulate(size)(i => 
      sum(glcm(i, ::).t :* DenseVector.ones[Double](size)) * pow(i - meanI, 2)
    )))
    val stdJ = sqrt(sum(DenseVector.tabulate(size)(j => 
      sum(glcm(::, j) :* DenseVector.ones[Double](size)) * pow(j - meanJ, 2)
    )))
    
    for {
      i <- 0 until size
      j <- 0 until size
      if glcm(i, j) > 0
    } {
      val value = glcm(i, j)
      contrast += value * pow(i - j, 2)
      correlation += value * (i - meanI) * (j - meanJ) / (stdI * stdJ)
      energy += value * value
      homogeneity += value / (1 + pow(i - j, 2))
      entropy -= value * log(value)
    }
    
    GLCMFeatures(contrast, correlation, energy, homogeneity, entropy)
  }

  def computeLBP(image: Array[Array[Int]], radius: Int = 1): Array[Array[Int]] = {
    val result = Array.ofDim[Int](image.length - 2 * radius, image(0).length - 2 * radius)
    
    for {
      x <- radius until image.length - radius
      y <- radius until image(0).length - radius
    } {
      var lbpValue = 0
      val centerValue = image(x)(y)
      
      // Sample points in clockwise order
      val points = for {
        i <- 0 until 8
        angle = 2 * math.Pi * i / 8
        px = (x + radius * math.cos(angle)).round.toInt
        py = (y + radius * math.sin(angle)).round.toInt
      } yield image(px)(py)
      
      // Compute LBP code
      for (i <- points.indices) {
        if (points(i) >= centerValue) {
          lbpValue |= (1 << i)
        }
      }
      
      result(x - radius)(y - radius) = lbpValue
    }
    
    result
  }

  def computeUniformLBP(image: Array[Array[Int]], radius: Int = 1): Array[Array[Int]] = {
    val basicLBP = computeLBP(image, radius)
    val uniformPatterns = getUniformPatterns()
    
    Array.tabulate(basicLBP.length, basicLBP(0).length) { (x, y) =>
      uniformPatterns.getOrElse(basicLBP(x)(y), 59) // Non-uniform patterns map to 59
    }
  }

  def computeLBPHistogram(lbpImage: Array[Array[Int]], numBins: Int = 60): Array[Double] = {
    val histogram = Array.fill(numBins)(0.0)
    val totalPixels = lbpImage.length * lbpImage(0).length
    
    for {
      x <- lbpImage.indices
      y <- lbpImage(0).indices
    } {
      histogram(lbpImage(x)(y)) += 1
    }
    
    // Normalize
    histogram.map(_ / totalPixels)
  }

  private def getUniformPatterns(): Map[Int, Int] = {
    val patterns = for {
      i <- 0 until 256
      binary = String.format("%8s", i.toBinaryString).replace(' ', '0')
      transitions = countTransitions(binary + binary(0))
      if transitions <= 2
    } yield i
    
    patterns.zipWithIndex.toMap
  }

  private def countTransitions(binary: String): Int = {
    var count = 0
    for (i <- 0 until binary.length - 1) {
      if (binary(i) != binary(i + 1)) count += 1
    }
    count
  }

  def computeHaralickFeatures(
    image: Array[Array[Int]], 
    distances: Array[Int] = Array(1, 2, 3), 
    angles: Array[Double] = Array(0.0, math.Pi/4, math.Pi/2, 3*math.Pi/4)
  ): Array[GLCMFeatures] = {
    val features = for {
      d <- distances
      a <- angles
    } yield {
      val glcm = computeGLCM(image, d, a)
      extractGLCMFeatures(glcm)
    }
    features.toArray
  }

  def computeMultiscaleLBP(
    image: Array[Array[Int]], 
    radii: Array[Int] = Array(1, 2, 3)
  ): Array[Array[Double]] = {
    radii.map { radius =>
      val lbpImage = computeUniformLBP(image, radius)
      computeLBPHistogram(lbpImage)
    }
  }
} 