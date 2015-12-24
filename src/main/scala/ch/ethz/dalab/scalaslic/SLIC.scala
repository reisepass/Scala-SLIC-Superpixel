package ch.ethz.dalab.scalaslic

import java.util.concurrent.atomic.AtomicInteger
import scala.collection.mutable.HashMap
import scala.collection.mutable.HashSet
import scala.collection.mutable.ListBuffer
import scala.collection.mutable.Stack
import scala.collection.parallel.mutable.ParMap
import scala.util.Random
import breeze.linalg.Vector
import breeze.numerics.floor
import breeze.numerics.round
import breeze.numerics.sqrt



case class DatumCord[DataCont](x: Int, y: Int, z: Int, cont: DataCont)

class SLICgreySimple(myImage: Array[Array[Array[Int]]], myS: Int, myM: Double, myMaxItr: Int = 30, myMinSize: Int = (-1)) extends SLIC[(Int)](
  distFn = (a: Int, b: Int) => Math.sqrt(Math.pow(a - b, 2)),
  rollingAvgFn = (a: Int, b: Int, n: Int) => {
    val sum = a * n + b
    (sum / (n + 1)).asInstanceOf[Int]
  },
  image = myImage,
  S = myS,
  inM = myM,
  user_cluster_normalization = false,
  maxIterations = myMaxItr,
  in_minSupSize = myMinSize) {

}

/*
 *   [DataType]  : Is the type which is stored in the grid structure. For grayscale images this could be an Int, For RGB images this could be (Int, Int, Int)
 *   distFn      : Defines how to calculate the distance between two datapoints. For RGB this could be sqrt( (r1-r2)^2  (g1 - g2)^2 + (b1-b2)^2)) for example 
 *   rollingAvgFn: Computes a rolling average of a series of datapoints. The use "  val newAvg = rollingAvgFn(old.cont, image(vX)(vY)(vZ), oldCount)" The old value, a new value to average in and how many points where used for the current average.
 *   image       : Input Image data in a grid structure 
 *   S           : Superpixel size 
 *   M           : Compactness, See paper for details. It weights the importance of spatial distance from the cluster center versus value or color distance. This parameter needs to be tuned carefully. 
 *   maxIterations: The maximum number of iterations before the result is returned 
 *   minChangePerIter: The minimum amount of cluster center change between interactions which will not cause a halt 
 *   user_cluster_normalization: Instead of setting M, it is possible to way spatial distance versus value distance by normalizing the color distance 
 *   in_minBlobSize: If superpixels end having a smaller size than this value they will be consolidated into a larger adjacent superpixel
 */
class SLIC[DataType](distFn: (DataType, DataType) => Double,
                     rollingAvgFn: ((DataType, DataType, Int) => DataType),
                     image: Array[Array[Array[DataType]]], 
                     S: Int, inM: Double = Double.MinValue,
                     maxIterations: Int = 15, minChangePerIter: Double = 0.000001,
                     user_cluster_normalization: Boolean = true, in_minSupSize: Int = (-1)) {
  
  

  //Initialize object variables and call center initialization functions
  val M = if (inM == Double.MinValue) S else inM
  val invwt = 1.0 / ((S / M) * (S / M))
  def distFn_m(a: DataType, b: DataType): Double = {     distFn(a, b)   }
  val xDim = image.size
  val yDim = image(0).size
  val zDim = image(0)(0).size
  val xS = if (xDim < S) xDim else S
  val yS = if (yDim < S) yDim else S
  val zS = if (zDim < S) zDim else S
  val minSupSize = if (in_minSupSize >= 0) in_minSupSize else (xS * yS * zS) / 2
  val dx = Array(-1, 1, 0, 0, 0, 0)
  val dy = Array(0, 0, -1, 1, 0, 0)
  val dz = Array(0, 0, 0, 0, -1, 1)
  val initCenter = placeInitSuperCenters(S)
  val centerCoords = moveInitSuperCentersToLowContrast(initCenter, 3)
  val clusterAssign = Array.fill(xDim, yDim, zDim) { -4 }.par
  println(" clusterAssign size = " + clusterAssign.size + " firstEl:" + clusterAssign(0)(0)(0))
  val centers = centerCoords.map(a => DatumCord(a._1, a._2, a._3, image(a._1)(a._2)(a._3))).par
  val lastDistance = Array.fill(xDim, yDim, zDim) { Double.MaxValue }.par
  val t0 = System.currentTimeMillis()
  val clusterMaxColDist = Array.fill(centers.size) { Double.MinValue }
  val t1 = System.currentTimeMillis()
  println((t0 - t1) + " for clusterMax Inst")
  val clusterMaxSpaceDist = Array.fill(centers.size) { Double.MinValue }
  println((System.currentTimeMillis() - t1) + "line 75")

  
  
  //Initializes the centers to a uniform grid
  def placeInitSuperCenters(gridInterval: Int): Array[(Int, Int, Int)] = {

    val xStart = if (xDim <= gridInterval) floor(xDim / 2) else round(gridInterval / 2)
    val yStart = if (yDim <= gridInterval) floor(yDim / 2) else round(gridInterval / 2)
    val zStart = if (zDim <= gridInterval) floor(zDim / 2) else round(gridInterval / 2)

    val out = for {
      x_s <- xStart until xDim by gridInterval;
      y_s <- yStart until yDim by gridInterval;
      z_s <- zStart until zDim by gridInterval
    } yield {
      (x_s, y_s, z_s)
    }

    out.toArray
  }

  //To prevent superpixel center initializations from being on a clear edge we perterb the centers to the lowest 'discrete' gradiant in a small radius around the original initilasation. 
  def moveInitSuperCentersToLowContrast(centers: Array[(Int, Int, Int)], purterbSpace: Int = 3): Array[(Int, Int, Int)] = {

    centers.map((a) => {
      val x = a._1
      val y = a._2
      val z = a._3
      var maxScore = 0.0
      var bestMove = (0, 0, 0)
      for {
        dx <- (-purterbSpace) until purterbSpace;
        dy <- (-purterbSpace) until purterbSpace;
        dz <- (-purterbSpace) until purterbSpace
      } {
        if (boundCheck(x + dx, y + dy, z + dz, xDim, yDim, zDim)) {
          val myCol = image(dx + x)(dy + y)(dz + z)
          var difSum = 0.0

          if (dx + x + 1 < xDim)
            difSum += distFn_m(myCol, image(dx + x + 1)(dy + y)(dz + z))
          if (dx + x - 1 >= 0)
            difSum += distFn_m(myCol, image(dx + x - 1)(dy + y)(dz + z))
          if (dy + y + 1 < yDim)
            difSum += distFn_m(myCol, image(dx + x)(dy + y + 1)(dz + z))
          if (dy + y - 1 >= 0)
            difSum += distFn_m(myCol, image(dx + x)(dy + y - 1)(dz + z))
          if (dz + z + 1 < zDim)
            difSum += distFn_m(myCol, image(dx + x)(dy + y)(dz + z + 1))
          if (dz + z - 1 >= 0)
            difSum += distFn_m(myCol, image(dx + x)(dy + y)(dz + z - 1))

          if (difSum < maxScore) {
            maxScore = difSum
            bestMove = (dx, dy, dz)
          }
        }
      }
      (x + bestMove._1, y + bestMove._2, z + bestMove._3)
    })

  }

  //Compute a datapoint's distance from a cluster. The point is passed by value the center is referred to by index because they are stored in a ParArray 
  def clusterDist(point: DatumCord[DataType], centerID: Int): Double = {
    val center = centers(centerID)
    val d_c = distFn_m(center.cont, point.cont)
    val d_s = sqrt(Math.pow(point.x - center.x, 2) + Math.pow(point.y - center.y, 2) + Math.pow(point.z - center.z, 2)) * invwt
    if (clusterMaxColDist(centerID) < d_c)
      clusterMaxColDist(centerID) = d_c
    if (clusterMaxSpaceDist(centerID) < d_s)
      clusterMaxSpaceDist(centerID) = d_s

    if (user_cluster_normalization)
      sqrt(Math.pow(d_c / clusterMaxColDist(centerID), 2) + Math.pow(d_s / clusterMaxSpaceDist(centerID), 2))
    else
      d_c + d_s

  }

  def boundCheck(x: Int, y: Int, z: Int, xDim: Int = xDim, yDim: Int = yDim, zDim: Int = zDim): Boolean = {
    val out = x >= 0 & x < xDim & y >= 0 & y < yDim & z >= 0 & z < zDim
    return out
  }

  //Returns the superpixel assignment per pixel. This is the only function which needs to be called after after initializing the object. 
  def calcSuperPixels(): Array[Array[Array[Int]]] = {

    assert(xDim > 0 & yDim > 0 & zDim > 0)
    val tt0 = System.currentTimeMillis()

    var lastChange = Double.MaxValue
    var otherHalt = 0
    println(centers.size + " centers.size")
    while (otherHalt < maxIterations & lastChange > minChangePerIter) {

      //Per Round
      val tPR = System.currentTimeMillis()

      (0 until centers.size).toList.par.map { cIDX =>
        {

          val tC = System.currentTimeMillis()
          val center = centers(cIDX)
          for (
            vX <- center.x - 2 * S until center.x + 2 * S;
            vY <- center.y - 2 * S until center.y + 2 * S;
            vZ <- center.z - 2 * S until center.z + 2 * S
          ) {
            if (boundCheck(vX, vY, vZ, xDim, yDim, zDim)) {

              val curVox = DatumCord(vX, vY, vZ, image(vX)(vY)(vZ))
              val curD = clusterDist(curVox, cIDX)
              if (lastDistance(vX)(vY)(vZ) > curD) {
                lastDistance(vX)(vY)(vZ) = curD
                clusterAssign(vX)(vY)(vZ) = cIDX
              }
            }
          }

        }
      }
      println((System.currentTimeMillis() - tPR) + " one round of cluste assignment updates")

      val updateRollingAvg = ParMap[Int, DatumCord[DataType]]()
      val updateCalcCount = ParMap[Int, Int]()

      val tMore = System.currentTimeMillis()

      println("")
      println("How Many Unlabled Voxels per slice")
      (0 until xDim).toList.map { vX => //TODO parallelize this.  (ParMap seems to not synchronize the way i expected it... Possible solution by using GroupBy)  
        {

          (0 until yDim).toList.map { vY =>
            {
              (0 until zDim).toList.map { vZ =>
                {

                  val usedForC = clusterAssign(vX)(vY)(vZ)
                  if (usedForC >= 0) {
                    if (updateRollingAvg.contains(usedForC)) {
                      val old = updateRollingAvg(usedForC)
                      val oldCount = updateCalcCount.get(usedForC).get
                      val newAvg = rollingAvgFn(old.cont, image(vX)(vY)(vZ), oldCount)
                      updateRollingAvg.put(usedForC, DatumCord(vX + old.x, vY + old.y, vZ + old.z, newAvg))
                      updateCalcCount.put(usedForC, oldCount + 1)
                    } else {
                      updateRollingAvg.put(usedForC, DatumCord(vX, vY, vZ, image(vX)(vY)(vZ)))
                      updateCalcCount.put(usedForC, 1)
                    }

                  }
                }

              }
            }
          }

        }

      }

      val lastCenters = centers.toArray
      updateRollingAvg.keys.foreach { key =>
        {
          val old = updateRollingAvg(key)
          val count = updateCalcCount(key)
          centers(key) = DatumCord(round(old.x / count), round(old.y / count), round(old.z / count), updateRollingAvg(key).cont)
        }
      }

      println()
      println((System.currentTimeMillis() - tMore) + " Update Cluster Centers [Centers.size= " + centers.size + "]")

      var totalChange = 0.0
      val tTime = System.currentTimeMillis()
      var howManyCountZero = 0;
      for (cIDX <- 0 until centers.size) {
        if (updateCalcCount.contains(cIDX)) {
          totalChange += clusterDist(lastCenters(cIDX), cIDX)
        } else {
          howManyCountZero += 1
        }
      }

      println((System.currentTimeMillis() - tTime) + " Error Detection")

      println("Centers Which where not used " + howManyCountZero)
      totalChange = totalChange / centers.size
      lastChange = totalChange

      otherHalt += 1
      println("Total Center Movement/centers.size " + totalChange)

    }

    val tConI = System.currentTimeMillis()
    val localC = clusterAssign.toArray
    val assingAfterEnf_i = enforceConnectivity_I(localC, minSupSize)
    println()
    println("Enfroce Connectivity (I " + (System.currentTimeMillis() - tConI) + ") ")
    println()
    return assingAfterEnf_i

  }

  //Center Assignments are not directly forced to form a continuous superpixel. The superpixels are made continuous by this function which reassigned isolated "islands" to the best adjacent superpixel.
  def enforceConnectivity_I(oldLabels: Array[Array[Array[Int]]], minBlobSize: Int): Array[Array[Array[Int]]] = {

    val xDim = oldLabels.size
    val yDim = oldLabels(0).size
    val zDim = oldLabels(0)(0).size
    val dx = Array(-1, 1, 0, 0, 0, 0)
    val dy = Array(0, 0, -1, 1, 0, 0)
    val dz = Array(0, 0, 0, 0, -1, 1)
    val outLabels = Array.fill(xDim, yDim, zDim) { -1 }

    val workingBlob = Stack[(Int, Int, Int)]()

    var newLabel = 0 //The old center id's are not needed b/c we we only use this fn as post processing 
    for (oX <- 0 until xDim; oY <- 0 until yDim; oZ <- 0 until zDim) {
      if (outLabels(oX)(oY)(oZ) < 0) {
        outLabels(oX)(oY)(oZ) = newLabel
        workingBlob.push((oX, oY, oZ))
        var adjLabel = -1

        for (n <- 0 until dx.size) {
          val curX = workingBlob.top._1 + dx(n)
          val curY = workingBlob.top._2 + dy(n)
          val curZ = workingBlob.top._3 + dz(n)
          if (boundCheck(curX, curY, curZ, xDim, yDim, zDim)) {
            if (outLabels(curX)(curY)(curZ) >= 0)
              adjLabel = outLabels(curX)(curY)(curZ)
          }

        }

        val totalBlob = Stack[(Int, Int, Int)]()
        while (workingBlob.size > 0) {
          val cur = workingBlob.pop
          totalBlob.push(cur)
          for (n <- 0 until dx.size) {
            val curX = cur._1 + dx(n)
            val curY = cur._2 + dy(n)
            val curZ = cur._3 + dz(n)
            if (boundCheck(curX, curY, curZ, xDim, yDim, zDim)) {
              if (outLabels(curX)(curY)(curZ) < 0 & oldLabels(curX)(curY)(curZ) == oldLabels(oX)(oY)(oZ)) {
                workingBlob.push((curX, curY, curZ))
                outLabels(curX)(curY)(curZ) = newLabel

              }

            }

          }
        }
        if (totalBlob.size <= minBlobSize) {
          if (adjLabel == (-1)) //This can happen if the first pixel is surrounded by pixels labeled differently
            totalBlob.foreach(a => { outLabels(a._1)(a._2)(a._3) = newLabel })
          else
            totalBlob.foreach(a => { outLabels(a._1)(a._2)(a._3) = adjLabel })
          newLabel -= 1
        }
        newLabel += 1
      }

    }

    return outLabels

  }

  //Simply returns a square grid as the superpixel mask. Could be used as a baseline comparison. 
  def calcSimpleSquaresSupPix(): Array[Array[Array[Int]]] = {
    assert(xDim > 0 & yDim > 0 & zDim > 0)
    val tt0 = System.currentTimeMillis()

    val maxSuperX = xDim / S
    val maxSuperY = yDim / S
    val maxSuperZ = zDim / S
    val supID = new AtomicInteger(0)
    val out = Array.fill(xDim, yDim, zDim) { -1 }
    for (X <- 0 to maxSuperX; Y <- 0 to maxSuperY; Z <- 0 to maxSuperZ) {
      var hasInc = false
      var curSupPix = supID.get

      val boundX = if (X == maxSuperX & (S + S * X < xDim)) xDim else S + S * X //Incase we cant divide the space into equaly sized squares we just extend the last square to fill the space
      val boundY = if (Y == maxSuperY & (S + S * Y < yDim)) yDim else S + S * Y
      val boundZ = if (Z == maxSuperZ & (S + S * Z < zDim)) zDim else S + S * Z
      for (x <- X * S to boundX; y <- Y * S to boundY; z <- Z * S to boundZ) {
        if (boundCheck(x, y, z, xDim, yDim, zDim)) {
          if (!hasInc) {
            curSupPix = supID.getAndIncrement
            hasInc = true
          }
          out(x)(y)(z) = curSupPix
          clusterAssign(x)(y)(z) = curSupPix
        }
      }

    }
    out
  }

  //Finds all pairs of superpixles which share an edge (Iterative version).  Simply loops over all pixels looking for two adjacent pixels which have different superpixel assigments in the mask array. 
  def findEdges_trueIter(supPixelIdmask: Array[Array[Array[Int]]], numSuperPix: Int): Array[scala.collection.mutable.Set[Int]] = {
    val xDim = supPixelIdmask.length
    val yDim = supPixelIdmask(0).length
    val zDim = supPixelIdmask(0)(0).length

    val conn = Array.fill(numSuperPix) { scala.collection.mutable.Set[Int]() }
    for (x <- 0 until xDim; y <- 0 until yDim; z <- 0 until zDim) {
      val me = supPixelIdmask(x)(y)(z)
      for {
        dx <- (-1) to 1;
        dy <- (-1) to 1;
        dz <- (-1) to 1
      } {
        if (boundCheck(x + dx, y + dy, z + dz, xDim, yDim, zDim)) {
          val other = supPixelIdmask(x + dx)(y + dy)(z + dz)

          if (me != other)
            conn(me) += (other)

        }
      }

    }
    conn
  }

  def findSupPixelBounds(supPixelId: Array[Array[Array[Int]]]): (Map[Int, List[DatumCord[DataType]]], Map[Int, HashMap[Int, Int]]) = {

    val cordWiseBlobs = HashMap[Int, List[DatumCord[DataType]]]()
    val supPixEdgeCount = HashMap[Int, HashMap[Int, Int]]()
    val xDim = supPixelId.length
    val yDim = supPixelId(0).length
    val zDim = supPixelId(0)(0).length

    for (x <- 0 until xDim; y <- 0 until yDim; z <- 0 until zDim) {
      val curId = supPixelId(x)(y)(z)
      val (allGood, myBlob, myEdges) = findBlobBounds_Rec(supPixelId, x, y, z, curId, new HashSet[(Int, Int, Int)](), new HashMap[Int, Int]())
      val blobWithData = myBlob.toList.map(a => DatumCord(a._1, a._2, a._3, image(a._1)(a._2)(a._3)))
      cordWiseBlobs.put(curId, blobWithData)
      supPixEdgeCount.put(curId, myEdges)
      if (curId == -1)
        print("should not find id = -1")
    }
    return (cordWiseBlobs.toMap, supPixEdgeCount.toMap)
  }

  //Finds all pairs of superpixles which share an edge (Recurssive Version). Performas an exaustive seearch for all superpixels which touch eachother starting from the inside and recursivly going out. 
  def findBlobBounds_Rec(supPixelId: Array[Array[Array[Int]]], x: Int, y: Int, z: Int, lastLabel: Int, myBlob: HashSet[(Int, Int, Int)], edgeCount: HashMap[Int, Int]): (Boolean, HashSet[(Int, Int, Int)], HashMap[Int, Int]) = {
    val xDim = supPixelId.length
    val yDim = supPixelId(0).length
    val zDim = supPixelId(0)(0).length

    if (x >= 0 & x < xDim & y >= 0 & y < yDim & z >= 0 & z < zDim) {

      if (supPixelId(x)(y)(z) == lastLabel) {
        myBlob += ((x, y, z))

        var outBlob = myBlob
        var outEdge = edgeCount
        for (d <- 0 until dx.size) {
          val nX = x + dx(d)
          val nY = y + dy(d)
          val nZ = z + dz(d)

          if (boundCheck(nX, nY, nZ, xDim, yDim, zDim))
            if (supPixelId(nX)(nY)(nZ) == lastLabel) {
              if (!(outBlob.contains((nX, nY, nZ)))) {
                val ret = findBlobBounds_Rec(supPixelId, nX, nY, nZ, lastLabel, outBlob, outEdge)
                if (ret._1) { //this boolean just checks if something was added
                  outBlob = ret._2
                  outEdge = ret._3
                }
              }
            } else {
              val old = outEdge.getOrElse(supPixelId(nX)(nY)(nZ), 0)
              outEdge.put(supPixelId(nX)(nY)(nZ), old + 1)
            }
        }
        return (true, outBlob, outEdge)
      } else
        return (false, myBlob, edgeCount)

    } else
      return (false, myBlob, edgeCount)

  }

  //Convert the superpixel mask of type Array[Array[Array[Int]]]  to something more easily usable for graph algorithms 
  def prepareGraph(supPixelId: Array[Array[Array[Int]]], featureFn: (List[DatumCord[DataType]]) => Vector[Double]): (List[Vector[Double]], HashMap[Int, Set[Int]]) = {
    import SLICutils.findEdges_simple;
    import SLICutils.findSupPixelCenterOfMassAndSize;
    val (supPixCenter, supPixSize) = findSupPixelCenterOfMassAndSize(supPixelId)
    val edgeMap = findEdges_simple(supPixelId, supPixCenter)
    val (cordWiseBlobs, supPixEdgeCount) = findSupPixelBounds(supPixelId)
    //make sure superPixelId's are ascending Ints
    val keys = cordWiseBlobs.keySet.toList.sorted
    val keymap = (keys.zip(0 until keys.size)).toMap
    def k(a: Int): Int = { keymap.get(a).get }

    val listOfFeatures = for (i <- 0 until keys.size) yield { featureFn(cordWiseBlobs.get(keys(i)).get) }
    val outEdgeMap = edgeMap.map(a => {
      val oldId = a._1
      val oldEdges = a._2
      (k(oldId), oldEdges.map { oldEdge => k(oldEdge) })
    })
    return (listOfFeatures.toList, outEdgeMap)

  }

  def clusterAssign(x: Int, y: Int, z: Int): Int = clusterAssign(x)(y)(z)

}

object SLICutils {
  val distFnCol = (a: (Int, Int, Int), b: (Int, Int, Int)) => sqrt(Math.pow(a._1 - b._1, 2) + Math.pow(a._2 - b._2, 2) + Math.pow(a._3 - b._3, 2))
  val sumFnCol = (a: (Int, Int, Int), b: (Int, Int, Int)) => ((a._1 + b._1, a._2 + b._2, a._3 + a._3))
  val normFnCol = (a: (Int, Int, Int), n: Int) => { ((a._1 / n, a._2 / n, a._3 / n)) }

  //Computes the size and center of mass for each superpixel 
  def findSupPixelCenterOfMassAndSize(supPixelId: Array[Array[Array[Int]]]): (HashMap[Int, (Int, Int, Int)], HashMap[Int, Int]) = {
    val xDim = supPixelId.length
    val yDim = supPixelId(0).length
    val zDim = supPixelId(0)(0).length
    val supPixCount = HashMap[Int, Int]()
    val supPixSum = HashMap[Int, (Int, Int, Int)]()
    for (x <- 0 until xDim; y <- 0 until yDim; z <- 0 until zDim) {
      val curId = supPixelId(x)(y)(z)
      val oldC = supPixCount.getOrElse(curId, 0)
      supPixCount.put(curId, oldC + 1)
      val oldS = supPixSum.getOrElse(curId, (0, 0, 0))
      supPixSum.put(curId, (oldS._1 + x, oldS._2 + y, oldS._3 + z))
    }
    supPixCount.keySet.foreach { key =>
      {
        val sum = supPixSum.get(key).get
        val cou = supPixCount.get(key).get
        supPixSum.put(key, (sum._1 / cou, sum._2 / cou, sum._3 / cou))
      }
    }
    return (supPixSum, supPixCount)
  }

  def findEdges_simple_array(supPixelIdmask: Array[Array[Array[Int]]], numSuperPix: Int): Array[scala.collection.mutable.Set[Int]] = {
    val (centers, size) = findSupPixelCenterOfMassAndSize(supPixelIdmask)
    val outAsMap = findEdges_simple(supPixelIdmask, centers)
    val numSup = outAsMap.keySet.size
    val keys = outAsMap.keySet.toArray.sorted

    val conn = Array.fill(numSuperPix) { scala.collection.mutable.Set[Int]() }
    keys.map(k => conn(k) = scala.collection.mutable.Set(outAsMap.get(k).get.toArray: _*))
    conn
  }

  //Quicky find edges between superpixels by moving in a straight line from the center in all directions until another superpixel is found. This will not find all edges or how large the edge is but it is very fast and in some cases enough.
  def findEdges_simple(supPixelIdmask: Array[Array[Array[Int]]], centers: HashMap[Int, (Int, Int, Int)]): HashMap[Int, Set[Int]] = {
    val xDim = supPixelIdmask.length
    val yDim = supPixelIdmask(0).length
    val zDim = supPixelIdmask(0)(0).length
    val out = HashMap[Int, Set[Int]]()
    centers.keySet.foreach { id =>
      {
        val c = centers.get(id).get
        val edges = ListBuffer[Int]()

        //Walk in each straight direction and find the next superPixel
        var dx = 1
        var stop = false
        while (dx + c._1 < xDim & !stop) {
          if (supPixelIdmask(c._1 + dx)(c._2)(c._3) != id) {
            edges += (supPixelIdmask(c._1 + dx)(c._2)(c._3))
            stop = true
          }
          dx += 1
        }
        stop = false
        dx = (-1)
        while (dx + c._1 >= 0 & !stop) {
          if (supPixelIdmask(c._1 + dx)(c._2)(c._3) != id) {
            edges += (supPixelIdmask(c._1 + dx)(c._2)(c._3))
            stop = true
          }
          dx -= 1
        }
        stop = false
        var dy = 1
        while (dy + c._2 < yDim & !stop) {
          if (supPixelIdmask(c._1)(c._2 + dy)(c._3) != id) {
            edges += (supPixelIdmask(c._1)(c._2 + dy)(c._3))
            stop = true
          }
          dy += 1
        }
        stop = false
        dy = (-1)
        while (dy + c._2 >= 0 & !stop) {
          if (supPixelIdmask(c._1)(c._2 + dy)(c._3) != id) {
            edges += (supPixelIdmask(c._1)(c._2 + dy)(c._3))
            stop = true
          }
          dy -= 1
        }
        stop = false
        var dz = 1
        while (dz + c._3 < zDim & !stop) {
          if (supPixelIdmask(c._1)(c._2)(c._3 + dz) != id) {
            edges += (supPixelIdmask(c._1)(c._2)(c._3 + dz))
            stop = true
          }
          dz += 1
        }
        stop = false
        dz = (-1)
        while (dz + c._3 >= 0 & !stop) {
          if (supPixelIdmask(c._1)(c._2)(c._3 + dz) != id) {
            edges += (supPixelIdmask(c._1)(c._2)(c._3 + dz))
            stop = true
          }
          dz -= 1
        }
        out.put(id, edges.toSet)

      }
    }
    return out
  }

}
