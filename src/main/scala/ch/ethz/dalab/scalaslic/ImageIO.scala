package ch.ethz.dalab.scalaslic

import java.awt.image.BufferedImage
import java.io.File
import javax.imageio.ImageIO
import scala.util.Try

object ImageIO {
  def readImage(path: String): Try[Array[Array[Array[Int]]]] = Try {
    val image = ImageIO.read(new File(path))
    val width = image.getWidth
    val height = image.getHeight
    val result = Array.ofDim[Int](width, height, 3)
    
    for {
      x <- 0 until width
      y <- 0 until height
    } {
      val rgb = image.getRGB(x, y)
      result(x)(y)(0) = (rgb >> 16) & 0xFF  // Red
      result(x)(y)(1) = (rgb >> 8) & 0xFF   // Green
      result(x)(y)(2) = rgb & 0xFF          // Blue
    }
    
    result
  }

  def writeImage(image: Array[Array[Array[Int]]], path: String): Try[Unit] = Try {
    val width = image.length
    val height = image(0).length
    val bufferedImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB)
    
    for {
      x <- 0 until width
      y <- 0 until height
    } {
      val r = image(x)(y)(0)
      val g = image(x)(y)(1)
      val b = image(x)(y)(2)
      val rgb = (r << 16) | (g << 8) | b
      bufferedImage.setRGB(x, y, rgb)
    }
    
    val format = path.substring(path.lastIndexOf('.') + 1)
    ImageIO.write(bufferedImage, format, new File(path))
  }

  def readGrayscaleImage(path: String): Try[Array[Array[Int]]] = Try {
    val image = ImageIO.read(new File(path))
    val width = image.getWidth
    val height = image.getHeight
    val result = Array.ofDim[Int](width, height)
    
    for {
      x <- 0 until width
      y <- 0 until height
    } {
      val rgb = image.getRGB(x, y)
      val r = (rgb >> 16) & 0xFF
      val g = (rgb >> 8) & 0xFF
      val b = rgb & 0xFF
      result(x)(y) = ((r * 0.299 + g * 0.587 + b * 0.114).toInt)
    }
    
    result
  }

  def writeGrayscaleImage(image: Array[Array[Int]], path: String): Try[Unit] = Try {
    val width = image.length
    val height = image(0).length
    val bufferedImage = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY)
    
    for {
      x <- 0 until width
      y <- 0 until height
    } {
      val gray = image(x)(y)
      val rgb = (gray << 16) | (gray << 8) | gray
      bufferedImage.setRGB(x, y, rgb)
    }
    
    val format = path.substring(path.lastIndexOf('.') + 1)
    ImageIO.write(bufferedImage, format, new File(path))
  }

  def visualizeSuperpixels(
    image: Array[Array[Array[Int]]], 
    superpixelMask: Array[Array[Int]], 
    boundaryColor: (Int, Int, Int) = (255, 0, 0)
  ): Array[Array[Array[Int]]] = {
    val result = Array.tabulate(image.length, image(0).length, image(0)(0).length) { (x, y, c) =>
      image(x)(y)(c)
    }
    
    // Draw boundaries
    for {
      x <- result.indices
      y <- result(0).indices
      if isBoundaryPixel(superpixelMask, x, y)
    } {
      result(x)(y)(0) = boundaryColor._1
      result(x)(y)(1) = boundaryColor._2
      result(x)(y)(2) = boundaryColor._3
    }
    
    result
  }

  def visualizeSegmentation(
    superpixelMask: Array[Array[Int]], 
    colorMap: Map[Int, (Int, Int, Int)] = null
  ): Array[Array[Array[Int]]] = {
    val labels = superpixelMask.flatten.distinct.sorted
    val colors = if (colorMap == null) {
      // Generate random colors for each label
      labels.map { label =>
        val random = new scala.util.Random(label)
        label -> (
          random.nextInt(256),
          random.nextInt(256),
          random.nextInt(256)
        )
      }.toMap
    } else colorMap
    
    Array.tabulate(superpixelMask.length, superpixelMask(0).length, 3) { (x, y, c) =>
      val label = superpixelMask(x)(y)
      val color = colors(label)
      c match {
        case 0 => color._1
        case 1 => color._2
        case 2 => color._3
      }
    }
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