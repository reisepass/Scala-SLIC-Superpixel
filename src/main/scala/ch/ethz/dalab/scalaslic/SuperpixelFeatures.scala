package ch.ethz.dalab.scalaslic

import breeze.linalg._
import breeze.numerics._
import breeze.stats._

case class SuperpixelFeatures(
  meanColor: Array[Double],
  colorHistogram: Array[Array[Double]],
  textureGradients: Array[Double],
  shape: ShapeDescriptor,
  centroid: (Double, Double),
  area: Int
)

case class ShapeDescriptor(
  compactness: Double,
  eccentricity: Double,
  orientation: Double,
  boundaryLength: Int
)

object SuperpixelFeatures {
  def extractFeatures(
    image: Array[Array[Array[Int]]], 
    superpixelMask: Array[Array[Int]], 
    superpixelId: Int,
    histogramBins: Int = 8
  ): SuperpixelFeatures = {
    
    // Collect all pixels belonging to the superpixel
    val pixels = for {
      x <- superpixelMask.indices
      y <- superpixelMask(0).indices
      if superpixelMask(x)(y) == superpixelId
    } yield (x, y)
    
    val area = pixels.length
    
    // Calculate centroid
    val centroid = pixels.foldLeft((0.0, 0.0)) { case ((sumX, sumY), (x, y)) =>
      (sumX + x, sumY + y)
    } match {
      case (sumX, sumY) => (sumX / area, sumY / area)
    }
    
    // Calculate mean color
    val meanColor = Array.fill(image(0)(0).length)(0.0)
    pixels.foreach { case (x, y) =>
      for (c <- meanColor.indices) {
        meanColor(c) += image(x)(y)(c)
      }
    }
    for (c <- meanColor.indices) {
      meanColor(c) /= area
    }
    
    // Calculate color histograms for each channel
    val colorHistogram = Array.fill(image(0)(0).length)(Array.fill(histogramBins)(0.0))
    pixels.foreach { case (x, y) =>
      for (c <- colorHistogram.indices) {
        val bin = (image(x)(y)(c) * histogramBins / 256).min(histogramBins - 1)
        colorHistogram(c)(bin) += 1
      }
    }
    // Normalize histograms
    for (c <- colorHistogram.indices) {
      val sum = colorHistogram(c).sum
      for (b <- colorHistogram(c).indices) {
        colorHistogram(c)(b) /= sum
      }
    }
    
    // Calculate texture gradients using Sobel
    val (gradX, gradY, magnitude) = EdgeDetector.sobelOperator(image)
    val textureGradients = pixels.map { case (x, y) => magnitude(x)(y) }.toArray
    
    // Calculate shape descriptors
    val shape = calculateShapeDescriptors(pixels.toArray, centroid, superpixelMask)
    
    SuperpixelFeatures(
      meanColor = meanColor,
      colorHistogram = colorHistogram,
      textureGradients = textureGradients,
      shape = shape,
      centroid = centroid,
      area = area
    )
  }
  
  private def calculateShapeDescriptors(
    pixels: Array[(Int, Int)], 
    centroid: (Double, Double),
    superpixelMask: Array[Array[Int]]
  ): ShapeDescriptor = {
    // Calculate covariance matrix for orientation and eccentricity
    val covMatrix = DenseMatrix.zeros[Double](2, 2)
    pixels.foreach { case (x, y) =>
      val dx = x - centroid._1
      val dy = y - centroid._2
      covMatrix(0, 0) += dx * dx
      covMatrix(0, 1) += dx * dy
      covMatrix(1, 0) += dx * dy
      covMatrix(1, 1) += dy * dy
    }
    covMatrix :/= pixels.length.toDouble
    
    // Calculate eigenvalues for eccentricity
    val eigs = eigSym(covMatrix)
    val eccentricity = sqrt(1.0 - eigs.eigenvalues(0) / eigs.eigenvalues(1))
    
    // Calculate orientation from eigenvectors
    val orientation = atan2(eigs.eigenvectors(1, 0), eigs.eigenvectors(0, 0))
    
    // Calculate compactness (area / perimeter^2)
    val boundaryPixels = pixels.count { case (x, y) =>
      val neighbors = for {
        dx <- -1 to 1
        dy <- -1 to 1
        if dx != 0 || dy != 0
        nx = x + dx
        ny = y + dy
        if nx >= 0 && nx < superpixelMask.length && ny >= 0 && ny < superpixelMask(0).length
      } yield superpixelMask(nx)(ny) != superpixelMask(x)(y)
      
      neighbors.exists(identity)
    }
    
    val compactness = 4.0 * math.Pi * pixels.length / (boundaryPixels * boundaryPixels)
    
    ShapeDescriptor(
      compactness = compactness,
      eccentricity = eccentricity,
      orientation = orientation,
      boundaryLength = boundaryPixels
    )
  }
  
  def computeSuperpixelSimilarity(f1: SuperpixelFeatures, f2: SuperpixelFeatures): Double = {
    // Color similarity (using mean color)
    val colorDiff = (f1.meanColor, f2.meanColor).zipped.map((a, b) => (a - b) * (a - b)).sum
    val colorSim = exp(-colorDiff / 100.0)
    
    // Histogram similarity (using chi-square distance)
    val histSim = (f1.colorHistogram, f2.colorHistogram).zipped.map { (h1, h2) =>
      val chiSquare = (h1, h2).zipped.map { (a, b) =>
        if (a + b > 0) (a - b) * (a - b) / (a + b) else 0.0
      }.sum
      exp(-chiSquare / 2.0)
    }.sum / f1.colorHistogram.length
    
    // Texture similarity
    val textureDiff = (f1.textureGradients, f2.textureGradients).zipped.map((a, b) => (a - b) * (a - b)).sum
    val textureSim = exp(-textureDiff / 1000.0)
    
    // Shape similarity
    val shapeSim = exp(-(
      math.pow(f1.shape.compactness - f2.shape.compactness, 2) +
      math.pow(f1.shape.eccentricity - f2.shape.eccentricity, 2) +
      math.pow(math.sin(f1.shape.orientation - f2.shape.orientation), 2)
    ) / 3.0)
    
    // Spatial proximity
    val spatialDiff = math.pow(f1.centroid._1 - f2.centroid._1, 2) + math.pow(f1.centroid._2 - f2.centroid._2, 2)
    val spatialSim = exp(-spatialDiff / (100.0 * 100.0))
    
    // Weighted combination
    0.3 * colorSim + 0.2 * histSim + 0.2 * textureSim + 0.15 * shapeSim + 0.15 * spatialSim
  }
} 