package ch.ethz.dalab.scalaslic

import breeze.linalg._
import breeze.numerics._
import scala.collection.mutable

object ImageMetrics {
  case class SegmentationMetrics(
    boundaryRecall: Double,
    boundaryPrecision: Double,
    fMeasure: Double,
    underSegmentation: Double,
    achievedCompactness: Double
  )

  def evaluateSegmentation(
    groundTruth: Array[Array[Int]], 
    segmentation: Array[Array[Int]], 
    tolerance: Int = 2
  ): SegmentationMetrics = {
    // Compute boundary maps
    val gtBoundaries = extractBoundaries(groundTruth)
    val segBoundaries = extractBoundaries(segmentation)
    
    // Compute boundary precision and recall
    val (recall, precision) = computeBoundaryPR(gtBoundaries, segBoundaries, tolerance)
    val fMeasure = if (recall + precision > 0) 2 * recall * precision / (recall + precision) else 0.0
    
    // Compute undersegmentation error
    val underSegmentation = computeUnderSegmentation(groundTruth, segmentation)
    
    // Compute achieved compactness
    val compactness = computeCompactness(segmentation)
    
    SegmentationMetrics(
      boundaryRecall = recall,
      boundaryPrecision = precision,
      fMeasure = fMeasure,
      underSegmentation = underSegmentation,
      achievedCompactness = compactness
    )
  }

  def computePSNR(original: Array[Array[Array[Int]]], processed: Array[Array[Array[Int]]]): Double = {
    require(original.length == processed.length && 
            original(0).length == processed(0).length && 
            original(0)(0).length == processed(0)(0).length)
    
    val mse = {
      var sum = 0.0
      var count = 0
      for {
        x <- original.indices
        y <- original(0).indices
        c <- original(0)(0).indices
      } {
        val diff = original(x)(y)(c) - processed(x)(y)(c)
        sum += diff * diff
        count += 1
      }
      sum / count
    }
    
    if (mse > 0) {
      10 * math.log10(255 * 255 / mse)
    } else {
      Double.PositiveInfinity
    }
  }

  def computeSSIM(
    original: Array[Array[Array[Int]]], 
    processed: Array[Array[Array[Int]]], 
    windowSize: Int = 8
  ): Double = {
    val c1 = math.pow(0.01 * 255, 2)
    val c2 = math.pow(0.03 * 255, 2)
    
    var ssimSum = 0.0
    var windowCount = 0
    
    for {
      x <- 0 until original.length - windowSize + 1 by windowSize
      y <- 0 until original(0).length - windowSize + 1 by windowSize
    } {
      // Extract window
      val window1 = extractWindow(original, x, y, windowSize)
      val window2 = extractWindow(processed, x, y, windowSize)
      
      // Compute statistics
      val mean1 = computeMean(window1)
      val mean2 = computeMean(window2)
      val variance1 = computeVariance(window1, mean1)
      val variance2 = computeVariance(window2, mean2)
      val covariance = computeCovariance(window1, window2, mean1, mean2)
      
      // Compute SSIM for window
      val numerator = (2 * mean1 * mean2 + c1) * (2 * covariance + c2)
      val denominator = (mean1 * mean1 + mean2 * mean2 + c1) * (variance1 + variance2 + c2)
      val ssim = numerator / denominator
      
      ssimSum += ssim
      windowCount += 1
    }
    
    ssimSum / windowCount
  }

  def computeSuperpixelCompactness(superpixelMask: Array[Array[Int]]): Array[Double] = {
    val labels = superpixelMask.flatten.distinct
    labels.map { label =>
      val pixels = for {
        x <- superpixelMask.indices
        y <- superpixelMask(0).indices
        if superpixelMask(x)(y) == label
      } yield (x, y)
      
      if (pixels.nonEmpty) {
        // Calculate centroid
        val centroid = pixels.foldLeft((0.0, 0.0)) { case ((sumX, sumY), (x, y)) =>
          (sumX + x, sumY + y)
        } match {
          case (sumX, sumY) => (sumX / pixels.length, sumY / pixels.length)
        }
        
        // Calculate average distance from centroid
        val avgDistance = pixels.map { case (x, y) =>
          math.sqrt(math.pow(x - centroid._1, 2) + math.pow(y - centroid._2, 2))
        }.sum / pixels.length
        
        // Calculate perimeter
        val perimeter = pixels.count { case (x, y) =>
          isBoundaryPixel(superpixelMask, x, y)
        }
        
        // Compactness = 4π * area / perimeter^2
        if (perimeter > 0) {
          4 * math.Pi * pixels.length / (perimeter * perimeter)
        } else 1.0
      } else 0.0
    }
  }

  private def extractBoundaries(mask: Array[Array[Int]]): Array[Array[Boolean]] = {
    Array.tabulate(mask.length, mask(0).length) { (x, y) =>
      isBoundaryPixel(mask, x, y)
    }
  }

  private def computeBoundaryPR(
    groundTruth: Array[Array[Boolean]], 
    segmentation: Array[Array[Boolean]], 
    tolerance: Int
  ): (Double, Double) = {
    var truePositives = 0
    var falsePositives = 0
    var falseNegatives = 0
    
    for {
      x <- groundTruth.indices
      y <- groundTruth(0).indices
    } {
      if (groundTruth(x)(y)) {
        // Check if there's a matching boundary pixel in segmentation within tolerance
        val found = (-tolerance to tolerance).exists { dx =>
          (-tolerance to tolerance).exists { dy =>
            val nx = x + dx
            val ny = y + dy
            nx >= 0 && nx < segmentation.length && 
            ny >= 0 && ny < segmentation(0).length && 
            segmentation(nx)(ny)
          }
        }
        if (found) truePositives += 1
        else falseNegatives += 1
      }
      
      if (segmentation(x)(y)) {
        // Check if there's a matching boundary pixel in ground truth within tolerance
        val found = (-tolerance to tolerance).exists { dx =>
          (-tolerance to tolerance).exists { dy =>
            val nx = x + dx
            val ny = y + dy
            nx >= 0 && nx < groundTruth.length && 
            ny >= 0 && ny < groundTruth(0).length && 
            groundTruth(nx)(ny)
          }
        }
        if (!found) falsePositives += 1
      }
    }
    
    val recall = if (truePositives + falseNegatives > 0) {
      truePositives.toDouble / (truePositives + falseNegatives)
    } else 0.0
    
    val precision = if (truePositives + falsePositives > 0) {
      truePositives.toDouble / (truePositives + falsePositives)
    } else 0.0
    
    (recall, precision)
  }

  private def computeUnderSegmentation(
    groundTruth: Array[Array[Int]], 
    segmentation: Array[Array[Int]]
  ): Double = {
    val gtLabels = groundTruth.flatten.distinct
    val segLabels = segmentation.flatten.distinct
    
    var totalError = 0.0
    var totalArea = 0
    
    for (gtLabel <- gtLabels) {
      val gtRegion = mutable.Set[(Int, Int)]()
      val overlapping = mutable.Map[Int, mutable.Set[(Int, Int)]]()
      
      // Find pixels belonging to ground truth region and their corresponding segmentation labels
      for {
        x <- groundTruth.indices
        y <- groundTruth(0).indices
        if groundTruth(x)(y) == gtLabel
      } {
        gtRegion.add((x, y))
        val segLabel = segmentation(x)(y)
        overlapping.getOrElseUpdate(segLabel, mutable.Set()).add((x, y))
      }
      
      // Calculate error for this ground truth region
      val regionError = overlapping.values.map(_.size).sum - gtRegion.size
      totalError += regionError
      totalArea += gtRegion.size
    }
    
    if (totalArea > 0) totalError / totalArea else 0.0
  }

  private def computeCompactness(segmentation: Array[Array[Int]]): Double = {
    val compactness = computeSuperpixelCompactness(segmentation)
    compactness.sum / compactness.length
  }

  private def extractWindow(
    image: Array[Array[Array[Int]]], 
    startX: Int, 
    startY: Int, 
    size: Int
  ): Array[Double] = {
    val window = new Array[Double](size * size * image(0)(0).length)
    var idx = 0
    for {
      x <- startX until startX + size
      y <- startY until startY + size
      c <- 0 until image(0)(0).length
    } {
      window(idx) = image(x)(y)(c)
      idx += 1
    }
    window
  }

  private def computeMean(window: Array[Double]): Double = {
    window.sum / window.length
  }

  private def computeVariance(window: Array[Double], mean: Double): Double = {
    window.map(x => math.pow(x - mean, 2)).sum / window.length
  }

  private def computeCovariance(
    window1: Array[Double], 
    window2: Array[Double], 
    mean1: Double, 
    mean2: Double
  ): Double = {
    (window1, window2).zipped.map { (x, y) =>
      (x - mean1) * (y - mean2)
    }.sum / window1.length
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