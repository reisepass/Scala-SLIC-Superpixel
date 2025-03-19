package ch.ethz.dalab.scalaslic

object Morphology {
  def dilate(image: Array[Array[Boolean]], kernelSize: Int): Array[Array[Boolean]] = {
    val result = Array.ofDim[Boolean](image.length, image(0).length)
    val radius = kernelSize / 2
    
    for {
      x <- 0 until image.length
      y <- 0 until image(0).length
    } {
      var found = false
      for {
        kx <- -radius to radius if !found
        ky <- -radius to radius if !found
      } {
        val nx = x + kx
        val ny = y + ky
        if (nx >= 0 && nx < image.length && ny >= 0 && ny < image(0).length) {
          if (image(nx)(ny)) {
            found = true
            result(x)(y) = true
          }
        }
      }
    }
    result
  }

  def erode(image: Array[Array[Boolean]], kernelSize: Int): Array[Array[Boolean]] = {
    val result = Array.ofDim[Boolean](image.length, image(0).length)
    val radius = kernelSize / 2
    
    for {
      x <- 0 until image.length
      y <- 0 until image(0).length
    } {
      var allTrue = true
      for {
        kx <- -radius to radius if allTrue
        ky <- -radius to radius if allTrue
      } {
        val nx = x + kx
        val ny = y + ky
        if (nx >= 0 && nx < image.length && ny >= 0 && ny < image(0).length) {
          if (!image(nx)(ny)) {
            allTrue = false
          }
        }
      }
      result(x)(y) = allTrue
    }
    result
  }

  def open(image: Array[Array[Boolean]], kernelSize: Int): Array[Array[Boolean]] = {
    dilate(erode(image, kernelSize), kernelSize)
  }

  def close(image: Array[Array[Boolean]], kernelSize: Int): Array[Array[Boolean]] = {
    erode(dilate(image, kernelSize), kernelSize)
  }

  def skeletonize(image: Array[Array[Boolean]]): Array[Array[Boolean]] = {
    val result = Array.tabulate(image.length, image(0).length)((x, y) => image(x)(y))
    var changed = true
    
    while (changed) {
      changed = false
      val toRemove = scala.collection.mutable.Set[(Int, Int)]()
      
      for {
        x <- 1 until image.length - 1
        y <- 1 until image(0).length - 1
        if result(x)(y)
      } {
        // Count neighbors
        var neighbors = 0
        var transitions = 0
        val neighborhood = Array.ofDim[Boolean](3, 3)
        
        for {
          i <- -1 to 1
          j <- -1 to 1
        } {
          neighborhood(i+1)(j+1) = result(x+i)(y+j)
          if (result(x+i)(y+j)) neighbors += 1
        }
        
        // Count transitions from 0 to 1 in clockwise order
        for {
          i <- 0 until 8
        } {
          val p1 = getPixel(neighborhood, i)
          val p2 = getPixel(neighborhood, (i + 1) % 8)
          if (!p1 && p2) transitions += 1
        }
        
        // Zhang-Suen thinning conditions
        if (neighbors >= 2 && neighbors <= 6 && 
            transitions == 1 && 
            !isEndpoint(neighborhood)) {
          toRemove.add((x, y))
          changed = true
        }
      }
      
      // Remove marked pixels
      for ((x, y) <- toRemove) {
        result(x)(y) = false
      }
    }
    
    result
  }
  
  private def getPixel(neighborhood: Array[Array[Boolean]], index: Int): Boolean = {
    index match {
      case 0 => neighborhood(0)(1) // N
      case 1 => neighborhood(0)(2) // NE
      case 2 => neighborhood(1)(2) // E
      case 3 => neighborhood(2)(2) // SE
      case 4 => neighborhood(2)(1) // S
      case 5 => neighborhood(2)(0) // SW
      case 6 => neighborhood(1)(0) // W
      case 7 => neighborhood(0)(0) // NW
      case _ => false
    }
  }
  
  private def isEndpoint(neighborhood: Array[Array[Boolean]]): Boolean = {
    var count = 0
    for {
      i <- 0 until 3
      j <- 0 until 3
      if !(i == 1 && j == 1)
    } {
      if (neighborhood(i)(j)) count += 1
    }
    count == 1
  }
} 