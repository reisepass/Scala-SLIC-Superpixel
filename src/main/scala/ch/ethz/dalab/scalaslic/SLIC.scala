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
//import scala.actors.threadpool.AtomicInteger

/*
 *  DataType for simple greyscale images will just be Double. for RGB it could be Int or (Int,Int,Int)  for example 
 *  dataDistanceFn just needs to define the distance between two colors. For grey scale we use linear euclidian distnace
 *  S is the grid size for the initial placement of super pixel centers. 
 *  
 */

//case class DatumCord[DataType](x: Int, y: Int, z: Int, cont: DataType)

case class DatumCord[DataCont](x:Int,y:Int,z:Int,cont:DataCont)

class SLIC[DataType](distFn: (DataType, DataType) => Double,
                     rollingAvgFn: ((DataType, DataType, Int) => DataType),
                     normFn: ((DataType, Int) => DataType),
                     image: Array[Array[Array[DataType]]], S: Int, K:Int=5, maxIterations: Int = 15, minChangePerIter: Double = 0.000001,
                     connectivityOption: String = "Functional", debug: Boolean = true) {

  //TODO try out different tif image stacks and see what dim2 is 

  val xDim = image.size
  val yDim = image(0).size
  val zDim = image(0)(0).size

  val xS = if (xDim < S) xDim else S
  val yS = if (yDim < S) yDim else S
  val zS = if (zDim < S) zDim else S
  val minBlobSize = (xS * yS * zS) / 2 //TODO check if this is the right minSize
  //cpp code looks like its using this:  floor(round((xDim*yDim*zDim)/K)/4)   
  //i think floor(x/4) is the same as x >> 2
  val dx = Array(-1, 1, 0, 0, 0, 0)
  val dy = Array(0, 0, -1, 1, 0, 0)
  val dz = Array(0, 0, 0, 0, -1, 1)

  val initCenter = placeInitSuperCenters(S)
  val centerCoords = moveInitSuperCentersToLowContrast(initCenter, 3)
  val clusterAssign = Array.fill(xDim, yDim, zDim) { -4 }.par //TODO Can be changed to Range for less memory usage
  println(" clusterAssign size = " + clusterAssign.size +" firstEl:"+clusterAssign(0)(0)(0))
  val centers = centerCoords.map(a => DatumCord(a._1, a._2, a._3, image(a._1)(a._2)(a._3))).par
  val lastDistance = Array.fill(xDim, yDim, zDim) { Double.MaxValue }.par
  val t0 = System.currentTimeMillis()
  val clusterMaxColDist = Array.fill(centers.size) { Double.MinValue }
  val t1 = System.currentTimeMillis()
  println((t0 - t1) + " for clusterMax Inst")
  val clusterMaxSpaceDist = Array.fill(centers.size) { Double.MinValue }
  println((System.currentTimeMillis() - t1) + "line 75")

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

  def moveInitSuperCentersToLowContrast(centers: Array[(Int, Int, Int)], purterbSpace: Int = 3): Array[(Int, Int, Int)] = {

    centers.map((a) => {
      val x = a._1
      val y = a._2
      val z = a._3
      var maxScore = Double.MaxValue
      var bestMove = (0, 0, 0)
      for {
        dx <- (-purterbSpace) until purterbSpace;
        dy <- (-purterbSpace) until purterbSpace;
        dz <- (-purterbSpace) until purterbSpace
      } {
        if (dx + x > 0 & dy + y > 0 & dz + z > 0 & dx + x < xDim & dy + y < yDim & dz + z < zDim) {
          val myCol = image(dx + x)(dy + y)(dz + z)
          var difSum = 0.0
          //TODO I am not checking bounds here b/c .getVoxel retunrs zero and what we want the behavior to be on the edges is not clear

          difSum += distFn(myCol, image(dx + x + 1)(dy + y)(dz + z))
          difSum += distFn(myCol, image(dx + x - 1)(dy + y)(dz + z))
          difSum += distFn(myCol, image(dx + x)(dy + y + 1)(dz + z))
          difSum += distFn(myCol, image(dx + x)(dy + y - 1)(dz + z))
          difSum += distFn(myCol, image(dx + x)(dy + y)(dz + z + 1))
          difSum += distFn(myCol, image(dx + x)(dy + y)(dz + z - 1))
          difSum
          if (difSum < maxScore) {
            maxScore = difSum
            bestMove = (dx, dy, dz)
          }
        }
      }
      (x + bestMove._1, y + bestMove._2, z + bestMove._3)
    })

  }

  def clusterDist(point: DatumCord[DataType], centerID: Int): Double = {
    val center = centers(centerID)
    val d_c = distFn(center.cont, point.cont)
    val d_s = sqrt(Math.pow(point.x - center.x, 2) + Math.pow(point.y - center.y, 2) + Math.pow(point.z - center.z, 2))
    if (clusterMaxColDist(centerID) < d_c)
      clusterMaxColDist(centerID) = d_c
    if (clusterMaxSpaceDist(centerID) < d_s)
      clusterMaxSpaceDist(centerID) = d_s
    //TODO question: should i update the clusterMax all the time or just after each round ? 
    sqrt(Math.pow(d_c / clusterMaxColDist(centerID), 2) + Math.pow(d_s / clusterMaxSpaceDist(centerID), 2))

  }

  def boundCheck(x: Int, y: Int, z: Int, xDim: Int = xDim, yDim: Int = yDim, zDim: Int = zDim): Boolean = {
    val out = x >= 0 & x < xDim & y >= 0 & y < yDim & z >= 0 & z < zDim
    return out
  }

  //This differs from findConnection_f in that it does not care about the true cluster label. It just starts with its own label 
  def findConnection_f2(x: Int, y: Int, z: Int, ourLabel: Int, myBlob: HashSet[(Int, Int, Int)], edgeCount: Vector[Int]): (Boolean, HashSet[(Int, Int, Int)], Vector[Int]) = {
    if (x >= 0 & x < xDim & y >= 0 & y < yDim & z >= 0 & z < zDim) {
      if (clusterAssign(x)(y)(z) == ourLabel) {

        if (myBlob.contains((x, y, z)))
          return (false, myBlob, edgeCount)
        else {

          myBlob += ((x, y, z))

          var outBlob = myBlob
          for (d <- 0 until dx.size) {
            val nX = x + dx(d)
            val nY = y + dy(d)
            val nZ = z + dz(d)

            if (boundCheck(nX, nY, nZ))
              if (clusterAssign(nX)(nY)(nZ) == ourLabel) {
                val ret = findConnection_f2(nX, nY, nZ, ourLabel, outBlob, edgeCount)
                if (ret._1) { //this boolean just checks if something was added
                  outBlob = ret._2
                }
              }
          }
          return (true, outBlob, edgeCount)
        }
      } else
        return (false, myBlob, edgeCount)
    } else
      return (false, myBlob, edgeCount)
  }

  def enforceConnectivity_F(clusterAssign: Array[Array[Array[Int]]], minBlobSize: Int): Array[Array[Array[Int]]] = {
    val blobConnections = Array.fill(xDim, yDim, zDim) { -1 }.par //TODO Can be changed to Range for less memory usage 

    //Current blog inuse is -3 
    //-1 is not seen yet 
    // 
    val maxLabel = centers.size //Anything above this means its a new blob 
    var lastAddedBlob = centers.size

    var edgeLabCount = Array.fill(centers.size) { 0 }
    var curBlob = new HashSet[(Int, Int, Int)]
    //We know these are connected to the center b/c they are the center
    for (i <- 0 until centers.size) {
      val cur = centers(i)

      blobConnections(cur.x)(cur.y)(cur.z) = i
    }

    val BREAK_SEARCH_AFTER = 20

    var debug_countEdge = 0
    def findConnection_f(x: Int, y: Int, z: Int, lastLabel: Int, myBlob: HashSet[(Int, Int, Int)], edgeCount: Vector[Int]): (Boolean, HashSet[(Int, Int, Int)], Vector[Int]) = {
      if (x >= 0 & x < xDim & y >= 0 & y < yDim & z >= 0 & z < zDim) {
        if (clusterAssign(x)(y)(z) == lastLabel) {

          // if (myBlob.contains((x,y,z)))
          if (blobConnections(x)(y)(z) == (-3))
            return (false, myBlob, edgeCount)
          else if (blobConnections(x)(y)(z) == lastLabel) {
            myBlob += ((x, y, z))
            return (true, myBlob, edgeCount)
          } else {

            // 

            myBlob += ((x, y, z))
            blobConnections(x)(y)(z) = (-3) //TODO Now that i added this back in I can change myBlob back into a list. B/c addin is faster for lists vs HashSets
            var outEdgeCount = edgeCount
            val n1 = findConnection_f(x + 1, y, z, lastLabel, myBlob, outEdgeCount)
            outEdgeCount = n1._3
            if (n1._1) {
              blobConnections(x)(y)(z) = lastLabel //this is kinda redundant from but if i ever take it parallel this could improve speed
              return (true, n1._2, n1._3)
            }

            val n2 = findConnection_f(x - 1, y, z, lastLabel, myBlob, outEdgeCount)
            outEdgeCount = n2._3
            if (n2._1) {
              blobConnections(x)(y)(z) = lastLabel

              return (true, n2._2, n2._3)
            }

            val n3 = findConnection_f(x, y + 1, z, lastLabel, myBlob, outEdgeCount)
            outEdgeCount = n3._3
            if (n3._1) {
              blobConnections(x)(y)(z) = lastLabel

              return (true, n3._2, n3._3)
            }

            val n4 = findConnection_f(x, y - 1, z, lastLabel, myBlob, outEdgeCount)
            outEdgeCount = n4._3
            if (n4._1) {
              blobConnections(x)(y)(z) = lastLabel

              return (true, n4._2, n4._3)
            }

            val n5 = findConnection_f(x, y, z + 1, lastLabel, myBlob, outEdgeCount)
            outEdgeCount = n5._3
            if (n5._1) {
              blobConnections(x)(y)(z) = lastLabel
              return (true, n5._2, n5._3)
            }

            val n6 = findConnection_f(x, y, z - 1, lastLabel, myBlob, outEdgeCount)
            outEdgeCount = n6._3
            if (n6._1) {
              blobConnections(x)(y)(z) = lastLabel
              return (true, n6._2, n6._3)
            } else
              return (false, myBlob, edgeCount)

          }
        } else {
          if (blobConnections(x)(y)(z) >= 0) {

            edgeCount(blobConnections(x)(y)(z)) += 1
            // debug_countEdge+=1
            return (false, myBlob, edgeCount)
          } else
            return (false, myBlob, edgeCount)

        }
      } else {
        return (false, myBlob, edgeCount)
      }

    }
    val newCenters = new ListBuffer[(DatumCord[DataType], Int)]

    val islandSize = new ListBuffer[Int]

    val remainingSmallIsland: scala.collection.parallel.mutable.ParSet[HashSet[(Int, Int, Int)]] = scala.collection.parallel.mutable.ParSet()
    val remainingBigIslands: scala.collection.parallel.mutable.ParSet[HashSet[(Int, Int, Int)]] = scala.collection.parallel.mutable.ParSet()

    val atomicNextLabel = new AtomicInteger(centers.size)
    Random.shuffle((0 until xDim).toList).map { vX => //would be nice to start each thread in a different section of the image to cause less (-3) collisions but preparing the coordinates takes alot of time so its not worth it 
      {
        (0 until yDim).toList.map { vY =>
          {
            (0 until zDim).toList.map { vZ =>
              if (blobConnections(vX)(vY)(vZ) == (-1)) {

                //val foundGroup = findConnection(vX, vY, vZ, clusterAssign(vX)(vY)(vZ))

                val altReturn = if (connectivityOption == "Functional2") findConnection_f2(vX, vY, vZ, clusterAssign(vX)(vY)(vZ), new HashSet[(Int, Int, Int)], Vector(Array.fill(centers.size) { 0 })) else findConnection_f(vX, vY, vZ, clusterAssign(vX)(vY)(vZ), new HashSet[(Int, Int, Int)], Vector(Array.fill(centers.size) { 0 }))
                val foundGroup = if (altReturn._1) clusterAssign(vX)(vY)(vZ) else -1
                val curBlob = altReturn._2
                val edgeLabCount = altReturn._3.toArray
                if (connectivityOption == "Functional2") {
                  curBlob.foreach(cord => {
                    blobConnections(cord._1)(cord._2)(cord._3) = (-3)
                  })
                  if (curBlob.size > minBlobSize) {
                    val myNewLabel = atomicNextLabel.getAndIncrement()
                    curBlob.foreach(cord => {
                      clusterAssign(cord._1)(cord._2)(cord._3) = myNewLabel
                      blobConnections(cord._1)(cord._2)(cord._3) = (myNewLabel)
                    })

                  } else {

                    print("(B" + curBlob.size + ")")

                    var neighCount = new HashMap[Int, Int]()
                    var bestLab = -1
                    var bestCount = 0
                    curBlob.foreach((cord) => {
                      for (n <- 0 until dx.size) { //TODO after we compare to c++ move this into the while loop and count all edges
                        val curX = cord._1 + dx(n)
                        val curY = cord._2 + dy(n)
                        val curZ = cord._3 + dz(n)
                        if (boundCheck(curX, curY, curZ, xDim, yDim, zDim)) {
                          if (!curBlob.contains((curX, curY, curZ))) {
                            val other = clusterAssign(curX)(curY)(curZ)
                            val old = neighCount.getOrElse(other, 0)
                            neighCount.put(other, old + 1)
                            if (old + 1 > bestCount) {
                              bestCount = old + 1
                              bestLab = other
                            }

                          }

                        }
                      }

                    })
                    assert(bestLab != (-1))
                    curBlob.foreach((cord) => {
                      clusterAssign(cord._1)(cord._2)(cord._3) = bestLab
                    })

                  }
                } else if (foundGroup == (-1)) {
                  islandSize += (curBlob.size)
                  if (curBlob.size > minBlobSize) { //The cluster is big enough to warrent it creating a new center
                    remainingBigIslands.+=(curBlob)
                  } else {
                    remainingSmallIsland.+=(curBlob)
                  }
                } else if (foundGroup == clusterAssign(vX)(vY)(vZ)) {
                  curBlob.foreach(a => {
                    val x = a._1
                    val y = a._2
                    val z = a._3
                    blobConnections(x)(y)(z) = foundGroup //TODO this is unnecessary alreayd done inside findConnection_f
                  })
                  //   curBlob = new HashSet[(Int, Int, Int)]
                } else {

                  println("Error, should not get here #DASDWQ")
                  assert(false)
                }
              }
            }
          }
        }
      }
    }

    val tBig = System.currentTimeMillis()
    var newLabel = atomicNextLabel.get
    remainingBigIslands.foreach(blob => {
      blob.foreach(cord => {
        clusterAssign(cord._1)(cord._2)(cord._3) = newLabel
        blobConnections(cord._1)(cord._2)(cord._3) = newLabel
      })
      newLabel += 1
    })

    val listSmallIslands = remainingSmallIsland.flatMap(cordBlob => { cordBlob.toList }).toArray
    val remaining = scala.collection.mutable.HashSet[Int]() ++= (0 until listSmallIslands.size).toList

    while (remaining.size > 0) {
      println("#remaining(" + remaining.size + ")")
      val removeFromRemaining = Stack[Int]()
      remaining.foreach { cordIDX =>
        {
          val cord = listSmallIslands(cordIDX)
          for (n <- 0 until dx.size) { //TODO after we compare to c++ move this into the while loop and count all edges
            val curX = cord._1 + dx(n)
            val curY = cord._2 + dy(n)
            val curZ = cord._3 + dz(n)
            if (boundCheck(curX, curY, curZ, xDim, yDim, zDim)) {
              val neigh = blobConnections(curX)(curY)(curZ)
              if (neigh >= 0) {
                blobConnections(cord._1)(cord._2)(cord._3) = neigh
                clusterAssign(cord._1)(cord._2)(cord._3) = neigh
                removeFromRemaining.push(cordIDX)
              }
            }
          }
        }
      }
      removeFromRemaining.foreach { rmCord => remaining.remove(rmCord) }
    }
    /*
      remainingSmallIsland.foreach( curBlob=>{
        if(curBlob.size>0) {
                    
                    var neighCount = new HashMap[Int,Int]()
                    var bestLab = -1
                    var bestCount = 0
                    curBlob.foreach((cord) => {
                      for (n <- 0 until dx.size) { //TODO after we compare to c++ move this into the while loop and count all edges
                        val curX = cord._1 + dx(n)
                        val curY = cord._2 + dy(n)
                        val curZ = cord._3 + dz(n)
                        if (boundCheck(curX, curY, curZ, xDim, yDim, zDim)) {
                          if (!curBlob.contains((curX, curY, curZ))){
                            val other =clusterAssign(curX)(curY)(curZ) 
                              val old = neighCount.getOrElse(other, 0)
                              neighCount.put(other,old+1)
                              if(old+1>bestCount){
                                bestCount=old+1
                                bestLab = other
                              }
                                
                          }
                            
                        }
                      }

                    })
                    if(bestLab==(-1)){
                      println("Cord[s="+curBlob.size+"]:")
                      curBlob.foreach((cord) => {
                        print( "("+cord._1+","+cord._2+","+cord._3+")")
                      })
                    assert(bestLab!=(-1))
                    }
                    curBlob.foreach((cord) => {
                      clusterAssign(cord._1)(cord._2)(cord._3) = bestLab
                    })
                    
                    
                  }

      })
    */

    return clusterAssign
  }

  def enforceConnectivity_I(oldLabels: Array[Array[Array[Int]]], minBlobSize: Int): Array[Array[Array[Int]]] = {

    //The below is designed to mick related C++ code and hance uses alot of state 
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

        for (n <- 0 until dx.size) { //TODO after we compare to c++ move this into the while loop and count all edges
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
          if(adjLabel==(-1)) //This can happen if the first pixel is surrounded by pixels labeled differently
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

  def calcSimpleSquaresSupPix():Array[Array[Array[Int]]] = {
    assert(xDim > 0 & yDim > 0 & zDim > 0)
    val tt0 = System.currentTimeMillis()
    
    val maxSuperX = xDim/S
    val maxSuperY = yDim/S
    val maxSuperZ = zDim/S
    val supID = new AtomicInteger(0)
    val out = Array.fill(xDim,yDim,zDim){-1}
    for ( X <- 0 to maxSuperX; Y <- 0 to maxSuperY ; Z <- 0 to maxSuperZ){
      var hasInc =false
      var curSupPix = supID.get
      
      val boundX = if(X==maxSuperX& (S+S*X<xDim) ) xDim else S+S*X //Incase we cant divide the space into equaly sized squares we just extend the last square to fill the space
      val boundY = if(Y==maxSuperY& (S+S*Y<yDim)) yDim else S+S*Y
      val boundZ = if(Z==maxSuperZ& (S+S*Z<zDim)) zDim else S+S*Z
      for(x <- X*S to boundX ; y<- Y*S to boundY; z <- Z*S to boundZ){
        if(boundCheck(x, y, z, xDim, yDim, zDim)){
          if(!hasInc){
            curSupPix=supID.getAndIncrement
            hasInc=true
          }
                  out(x)(y)(z)=curSupPix
        clusterAssign(x)(y)(z)=curSupPix
        }
      }
        
    }
    out
  }
  
  def calcSuperPixels(): Array[Array[Array[Int]]] = {

    if (debug)
      println("(" + xDim + "," + yDim + "," + zDim + ") image Dims")
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
            vX <- center.x - 2 * S until center.x + 2 * S; //TODO change back to 2*S
            vY <- center.y - 2 * S until center.y + 2 * S;
            vZ <- center.z - 2 * S until center.z + 2 * S
          ) {
            if (boundCheck(vX,vY,vZ, xDim, yDim, zDim)) {

              val curVox = DatumCord(vX, vY, vZ, image(vX)(vY)(vZ))
              val curD = clusterDist(curVox, cIDX)
              if (lastDistance(vX)(vY)(vZ) > curD) {
                lastDistance(vX)(vY)(vZ) = curD
                clusterAssign(vX)(vY)(vZ) = cIDX
              }
            }
          }
          // println((System.currentTimeMillis()-tC)+" Just one Center "+cIDX) 
        }
      }
      println((System.currentTimeMillis() - tPR) + " one round of cluste assignment updates")

      val updateRollingAvg = ParMap[Int, DatumCord[DataType]]()
      val updateCalcCount = ParMap[Int, Int]()

      val tMore = System.currentTimeMillis()

      println("")
      println("How Many Unlabled Voxels per slice")
      (0 until xDim).toList.map { vX => //TODO parallize this.  (ParMap seems to not syncronize the way i expected it...)  
        { // Use groupBy  function then reduce to make it nice 
          var howManyNegOne = 0
          (0 until yDim).toList.map { vY =>
            {
              (0 until zDim).toList.map { vZ =>
                {

                  val usedForC = clusterAssign(vX)(vY)(vZ)
                  if (usedForC >=0) {
                    if (updateRollingAvg.contains(usedForC)) {
                      val old = updateRollingAvg(usedForC)
                      val oldCount = updateCalcCount.get(usedForC).get
                      val newAvg = rollingAvgFn( old.cont,image(vX)(vY)(vZ),oldCount)
                      updateRollingAvg.put(usedForC, DatumCord(vX + old.x, vY + old.y, vZ + old.z,newAvg ))
                      updateCalcCount.put(usedForC, oldCount + 1)
                    } else {
                      updateRollingAvg.put(usedForC, DatumCord(vX, vY, vZ, image(vX)(vY)(vZ)))
                      updateCalcCount.put(usedForC, 1)
                    }

                  } else {
                    howManyNegOne += 1
                    //(-1,null)
                  }
                }

              }
            }
          }
          print(" (" + howManyNegOne + ")")

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

    if (debug) {
      val tConI = System.currentTimeMillis()
      val localClust = clusterAssign.toArray
      val assingAfterEnf_i = enforceConnectivity_I(localClust, minBlobSize)
      val tConF = System.currentTimeMillis()
      println(" i fin ")
      val assignAfterEnf_f = enforceConnectivity_F(localClust, minBlobSize)
      println("Enfroce Connectivity (I " + (tConF - tConI) + ") (F " + (System.currentTimeMillis() - tConF) + ")")
      // val debugf2 = enforceConnectivity_F(assignAfterEnf_f,minBlobSize)
      return assignAfterEnf_f
    }
    if (connectivityOption == "Functional" || connectivityOption == "Functional2") {
      val tConF = System.currentTimeMillis()
      val localC = clusterAssign.toArray
      val assignAfterEnf_f = enforceConnectivity_F(localC, minBlobSize)
      println()
      println("Enfroce Connectivity (F " + (System.currentTimeMillis() - tConF) + ") ")
      println()
      return assignAfterEnf_f
    }
    if (connectivityOption == "Imperative") {
      val tConI = System.currentTimeMillis()
      val localC = clusterAssign.toArray
      val assingAfterEnf_i = enforceConnectivity_I(localC, minBlobSize)
      println()
      println("Enfroce Connectivity (I " + (System.currentTimeMillis() - tConI) + ") ")
      println()
      return assingAfterEnf_i

    }
    if (connectivityOption == "None") {
      return clusterAssign.toArray
    }

    return clusterAssign.toArray //TODO This should never be touched, I can force this by putting an assurt that the connectivityOption is only one of the above options
  }

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
      if(curId == -1)
        print("should not find id = -1")
    }
    return (cordWiseBlobs.toMap, supPixEdgeCount.toMap)
  }
  
  def findSupPixelBounds_I(supPixelId: Array[Array[Array[Int]]]): Map[Int, List[DatumCord[DataType]]] = {
        val xDim = supPixelId.length
    val yDim = supPixelId(0).length
    val zDim = supPixelId(0)(0).length
    val cordWiseBlobs = HashMap[Int, ListBuffer[DatumCord[DataType]]]()
    
    for (x <- 0 until xDim; y <- 0 until yDim; z <- 0 until zDim) {
      val curId = supPixelId(x)(y)(z)
            if(curId == -1)
        print("should not find id = -1")
        
      if(!cordWiseBlobs.contains(curId))
        cordWiseBlobs.put(curId,new ListBuffer[DatumCord[DataType]]())
      cordWiseBlobs.get(curId).get+=DatumCord(x,y,z,image(x)(y)(z))
    }
    
    cordWiseBlobs.toMap.map{ case (key,value)=>(key,value.toList)}
  }

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

  def prepareGraph(supPixelId: Array[Array[Array[Int]]], featureFn: (List[DatumCord[DataType]]) => Vector[Double]): (List[Vector[Double]], HashMap[Int, Set[Int]]) = {
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