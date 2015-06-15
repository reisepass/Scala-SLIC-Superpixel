import ij._

import ij.io.Opener
import javax.imageio.ImageIO
import java.io.File
import java.io.ByteArrayOutputStream
import breeze.linalg._
import breeze.stats.DescriptiveStats._
import breeze.stats._
import breeze.numerics._
import org.apache.commons.math3.ml.distance.DistanceMeasure
import ij.plugin.Duplicator
import scala.collection.mutable.HashSet
import scala.collection.mutable.ListBuffer
import scala.collection.immutable.HashMap
import scala.collection.mutable.LinkedList
import scala.util.Random
import ch.ethz.dalab.scalaslic.SLIC
import ch.ethz.dalab.scalaslic.DatumCord

object testImageJ {

  def main(args: Array[String]): Unit = {

    val options: Map[String, String] = args.map { arg =>
      arg.dropWhile(_ == '-').split('=') match {
        case Array(opt, v) => (opt -> v)
        case Array(opt)    => (opt -> "true")
        case _             => throw new IllegalArgumentException("Invalid argument: " + arg)
      }
    }.toMap

    val tft0 = System.currentTimeMillis()
    // val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/testing_groundtruth.tif"
    //val examplePath = "/home/mort/workspace/back_ScalaSLIC/data/training.tif"
    val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/training_80_84_165.tif"
    // val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/training_48_40_40.tif"
    //val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/training_9_9_10.tif"
    // val examplePath = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/training_saveAsimgJ.tif"
    val examplePa3 = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/testing_groundtruth_first5.tif"
    val examplePath2 = "/home/mort/workspace/Scala-SLIC-Superpixel/ScalaSLIC/data/training.tif"
    val S = 20 // Superpixel center grid interval 

    val opener = new Opener();
    val imp2 = opener.openImage(examplePath);

    val tft1 = System.currentTimeMillis()
    println((tft1 - tft0) + " OpenImage() time")
    //imp2.show() 

    val aStack = imp2.getStack

    val v1 = aStack.getVoxel(1, 1, 1)
    val v2 = aStack.getVoxel(1, 2, 1)
    val v3 = aStack.getVoxel(1, 1, 2)
    val bitDepa = aStack.getBitDepth()
    val colModa = aStack.getColorModel()
    val zsize = aStack.getSize //Zsize 
    val allTheD = imp2.getDimensions
    val xDim = allTheD(0)
    val yDim = allTheD(1)
    val zDim = allTheD(3)
    println("(" + xDim + "," + yDim + "," + zDim + ") image Dims")
    assert(xDim > 0 & yDim > 0 & zDim > 0)

    val tFu = System.currentTimeMillis()
    val copyImg = Array.fill(xDim,yDim,zDim)(0.0)
    for(x <- 0 until xDim; y <- 0 until yDim; z <- 0 until zDim){
      copyImg(x)(y)(z)=aStack.getVoxel(x,y,z)
    }
    println(" copying image over took "+(System.currentTimeMillis()-tFu))
    
    /*
     * 
     * Test Some More 
     * 
     * 
     */
    val distFn = (a:Double,b:Double)=>sqrt(Math.pow(a-b,2))
    val sumFn = (a:Double,b:Double)=>(a+b)
    val normFn = (a:Double,b:Int)=>(a/b)
    
     
     val allOp =  new SLIC[Double](distFn,sumFn,normFn,copyImg,20,15)
    
    
    
    
    /*
     * 
     * Testing stuff 
     */
     
    val img3 = opener.openImage("/home/mort/workspace/dissolve-struct/data/generated/MSRC_ObjCategImageDatabase_v2/Images/1_10_s.bmp")
    //val aldim = img3.getDimensions
    //val pixF = img3.getPixel(0,0)
    
    val stack3 = img3.getStack()
    val x2 = stack3.getWidth
    val y2 = stack3.getHeight
    val z2 = stack3.getSize
    val rgbja = stack3.isRGB 
    val bitDep = stack3.getBitDepth()
    val colMod = stack3.getColorModel()
    val v21 =   stack3.getVoxel(x2-1, y2-2, z2-1)
    val v22 =   stack3.getVoxel(x2, y2, z2)
    val v223 =  stack3.getVoxel(x2, y2,0)
    val v23 =   stack3.getVoxel(1, 1, 2)
    val vwha =  stack3.getVoxel(1,1,0)
    val vwha2 = stack3.getVoxel(1,1,4)
    val vwha3 = stack3.getVoxel(1,1,3)
    val vwha4 = stack3.getVoxel(1,1,59)
    val vwha5 = stack3.getVoxel(1,2,0)
    val vwha6 = stack3.getVoxel(2,1,0)
    val vwha62 = stack3.getVoxel(0,0,0)
   val grr = colMod.getGreen(vwha62.asInstanceOf[Int])
    println("fun times")
    
    
    
    
    case class Grey3dCord(x: Int, y: Int, z: Int, grey: Double)
    val test = 5

    def placeInitSuperCenters(gridInterval: Int): Array[(Int, Int, Int)] = { //THis function may be unnecessary
      val out = for {
        x_s <- round(gridInterval / 2) until xDim by gridInterval;
        y_s <- round(gridInterval / 2) until yDim by gridInterval;
        z_s <- round(gridInterval / 2) until zDim by gridInterval
      } yield {
        (x_s, y_s, z_s)
      }
      out.toArray
    }

    def moveInitSuperCentersToLowContrast(centers: Array[(Int, Int, Int)], myStack: ImageStack, purterbSpace: Int = 3): Array[(Int, Int, Int)] = { //THis function may be unnecessary

      val xDim = myStack.getWidth()
      val yDim = myStack.getHeight()
      val zDim = myStack.getSize
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
            val myCol = myStack.getVoxel(dx + x, dy + y, dz + z)
            var difSum = 0.0
            //TODO I am not checking bounds here b/c .getVoxel retunrs zero and what we want the behavior to be on the edges is not clear
            difSum += pow((myCol - myStack.getVoxel(dx + x + 1, dy + y, dz + z)), 2)
            difSum += pow((myCol - myStack.getVoxel(dx + x - 1, dy + y, dz + z)), 2)
            difSum += pow((myCol - myStack.getVoxel(dx + x, dy + y + 1, dz + z)), 2)
            difSum += pow((myCol - myStack.getVoxel(dx + x, dy + y - 1, dz + z)), 2)
            difSum += pow((myCol - myStack.getVoxel(dx + x, dy + y, dz + z + 1)), 2)
            difSum += pow((myCol - myStack.getVoxel(dx + x, dy + y, dz + z - 1)), 2) //TODO add diagonal edge differences too
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

    val tt0 = System.currentTimeMillis()

    val initCenter = placeInitSuperCenters(S)
    val centerCoords = moveInitSuperCentersToLowContrast(initCenter, aStack, 3)
    val clusterAssign = Array.fill(xDim, yDim, zDim) { -1 }.par //TODO Can be changed to Range for less memory usage
    println(" clusterAssign size = " + clusterAssign.size + clusterAssign.toString())
    val lastDistance = Array.fill(xDim, yDim, zDim) { Double.MaxValue }.par
    val centers = centerCoords.map(a => Grey3dCord(a._1, a._2, a._3, aStack.getVoxel(a._1, a._2, a._3))).par //This should be init to the "lowest gradient 2x2x2 neearbye"

    val t0 = System.currentTimeMillis()
    println((t0 - tt0) + " 58 to 66")
    val clusterMaxColDist = Array.fill(centers.size) { Double.MinValue }
    val t1 = System.currentTimeMillis()
    println((t0 - t1) + " for clusterMax Inst")
    val clusterMaxSpaceDist = Array.fill(centers.size) { Double.MinValue }
    println((System.currentTimeMillis() - t1) + "line 75")

    /*
     * 
     * QUICK EXP on the retrival of voxels 
     */
    val random = new scala.util.Random(0)
    val tstack = System.currentTimeMillis()
    val totalRet = 100000
    for (i <- 0 until totalRet) {
      val ass = aStack.getVoxel(random.nextInt(xDim), random.nextInt(yDim), random.nextInt(zDim))

    }
    println("GetVoxel" + totalRet + " times takes: " + (System.currentTimeMillis() - tstack))
    val tClust = System.currentTimeMillis()
    for (i <- 0 until totalRet) {
      val ass = clusterAssign(random.nextInt(xDim))(random.nextInt(yDim))(random.nextInt(zDim))
    }
    println("matrix retrival " + totalRet + " times takes: " + (System.currentTimeMillis() - tClust))

    def clusterDist(point: Grey3dCord, centerID: Int): Double = {
      val center = centers(centerID)
      val d_c = sqrt(Math.pow(point.grey - center.grey, 2))
      val d_s = sqrt(Math.pow(point.x - center.x, 2) + Math.pow(point.y - center.y, 2) + Math.pow(point.z - center.z, 2))
      if (clusterMaxColDist(centerID) < d_c)
        clusterMaxColDist(centerID) = d_c
      if (clusterMaxSpaceDist(centerID) < d_s)
        clusterMaxSpaceDist(centerID) = d_s
      //TODO question: should i update the clusterMax all the time or just after each round ? 
      sqrt(Math.pow(d_c / clusterMaxColDist(centerID), 2) + Math.pow(d_s / clusterMaxSpaceDist(centerID), 2))

    }

    def distanceMeasure(left: Grey3dCord, right: Grey3dCord, clusterMaxC: Double, clusterMaxS: Double): Double = {
      val d_c = sqrt(Math.pow(left.grey - right.grey, 2))
      val d_s = sqrt(Math.pow(left.x - right.x, 2) + Math.pow(left.y - right.y, 2) + Math.pow(left.z - right.z, 2))
      sqrt(Math.pow(d_c / clusterMaxC, 2) + Math.pow(d_s / clusterMaxS, 2))
    }

    def enforceConnectedness() {
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

      def averageGrey3dCord(in: List[Grey3dCord]): Grey3dCord = {

        val avgX = round(in.map { a => a.x }.sum / in.size)
        val avgY = round(in.map { a => a.y }.sum / in.size)
        val avgZ = round(in.map { a => a.z }.sum / in.size)
        val avgGrey = in.map { a => a.grey }.sum / in.size
        return Grey3dCord(avgX, avgY, avgZ, avgGrey)

      }

      val BREAK_SEARCH_AFTER = 20

      var debug_countEdge=0
      def findConnection_f(x: Int, y: Int, z: Int, lastLabel: Int, myBlob: HashSet[(Int, Int, Int)],edgeCount:Vector[Int]): (Boolean, HashSet[(Int, Int, Int)],Vector[Int]) = {
        if (x >= 0 & x < xDim & y >= 0 & y < yDim & z >= 0 & z < zDim) {
          if (clusterAssign(x)(y)(z) == lastLabel) {

           // if (myBlob.contains((x,y,z)))
             if (blobConnections(x)(y)(z) == (-3))
              return (false, myBlob,edgeCount)
            else if (blobConnections(x)(y)(z) == lastLabel) {
              myBlob += ((x, y, z))
              return (true, myBlob,edgeCount)
            } else {

              // 
              
              myBlob += ((x, y, z))
              blobConnections(x)(y)(z) = (-3)//TODO Now that i added this back in I can change myBlob back into a list. B/c addin is faster for lists vs HashSets
              var outEdgeCount=edgeCount
              val n1 = findConnection_f(x + 1, y, z, lastLabel, myBlob,outEdgeCount)
              outEdgeCount = n1._3
              if (n1._1) {
                blobConnections(x)(y)(z) = lastLabel  //this is kinda redundant from but if i ever take it parallel this could improve speed
                return (true, n1._2,n1._3)
              }

              val n2 = findConnection_f(x - 1, y, z, lastLabel, myBlob,outEdgeCount)
              outEdgeCount= n2._3
              if(n2._1){
                blobConnections(x)(y)(z) = lastLabel
              
                return (true, n2._2,n2._3)
              }

              val n3 = findConnection_f(x, y + 1, z, lastLabel, myBlob,outEdgeCount)
              outEdgeCount= n3._3
              if(n3._1){
                blobConnections(x)(y)(z) = lastLabel
              
                return (true, n3._2,n3._3)
              }

              val n4 = findConnection_f(x, y - 1, z, lastLabel, myBlob,outEdgeCount)
              outEdgeCount= n4._3
              if(n4._1){
                blobConnections(x)(y)(z) = lastLabel
              
                return (true, n4._2,n4._3)
              }

              val n5 = findConnection_f(x, y, z + 1, lastLabel, myBlob,outEdgeCount)
              outEdgeCount= n5._3
              if(n5._1){
                blobConnections(x)(y)(z) = lastLabel
                return (true, n5._2,n5._3)
              }

              val n6 = findConnection_f(x, y, z - 1, lastLabel, myBlob,outEdgeCount)
              outEdgeCount= n6._3
              if (n6._1){
                blobConnections(x)(y)(z) = lastLabel
                return (true, n6._2,n6._3)
              }
              else
                return (false, myBlob,edgeCount)

            }
          } else {
            if( blobConnections(x)(y)(z)>=0){
             
             edgeCount(blobConnections(x)(y)(z))+=1
            // debug_countEdge+=1
             return (false, myBlob,edgeCount)  
            }
            else
              return (false, myBlob,edgeCount)              
                                                      
          }
        } else {
          return (false, myBlob,edgeCount)
        }

      }

      def findConnection(x: Int, y: Int, z: Int, lastLabel: Int): Int = { //TODO this could be changed to a Boolean return found or not 
        if (x >= 0 & x < xDim & y >= 0 & y < yDim & z >= 0 & z < zDim) {
          if (clusterAssign(x)(y)(z) == lastLabel) {

            if (blobConnections(x)(y)(z) == (-3))
              return -1
            else if (blobConnections(x)(y)(z) == lastLabel) {

              return lastLabel
            } else {
              curBlob += ((x, y, z))
              blobConnections(x)(y)(z) = (-3)

              val n1 = findConnection(x + 1, y, z, lastLabel)
              
              if (n1 == lastLabel)
                return lastLabel

              val n2 = findConnection(x - 1, y, z, lastLabel)
              
              if (n2 == lastLabel)
                return lastLabel

              val n3 = findConnection(x, y + 1, z, lastLabel)
              
              if (n3 == lastLabel)
                return lastLabel

              val n4 = findConnection(x, y - 1, z, lastLabel)
              
              if (n4 == lastLabel)
                return lastLabel

              val n5 = findConnection(x, y, z + 1, lastLabel)
              
              if (n5 == lastLabel)
                return lastLabel

              val n6 = findConnection(x, y, z - 1, lastLabel)
              
              if (n6 == lastLabel)
                return lastLabel
              else
                return -1

            }
          } else {
            edgeLabCount(clusterAssign(x)(y)(z)) += 1
            return -1
          }
        } else {
          return -1
        }

      }

      val newCenters = new ListBuffer[(Grey3dCord, Int)]

      val islandSize = new ListBuffer[Int]
      
      val what = List(1,2,3,4).par
      val remainingIslandBlobs : scala.collection.parallel.mutable.ParSet[HashSet[(Int, Int, Int)]]= scala.collection.parallel.mutable.ParSet()
      
      
       Random.shuffle((0 until xDim).toList).par.map { vX => //would be nice to start each thread in a different section of the image to cause less (-3) collisions but preparing the coordinates takes alot of time so its not worth it 
        { 
          (0 until yDim).toList.map { vY =>
            {
              (0 until zDim).toList.map { vZ =>
                if (blobConnections(vX)(vY)(vZ) == (-1)) {
                  
                  //val foundGroup = findConnection(vX, vY, vZ, clusterAssign(vX)(vY)(vZ))
                  
                  val altReturn = findConnection_f(vX,vY,vZ,clusterAssign(vX)(vY)(vZ),new HashSet[(Int,Int,Int)],Vector(Array.fill(centers.size){0}))
                  val foundGroup = if(altReturn._1) clusterAssign(vX)(vY)(vZ) else -1
                  val curBlob = altReturn._2
                  val edgeLabCount = altReturn._3.toArray
                 // assert(edgeMap.keySet.size >0)//TODO you need to deal with this case.
                  if (foundGroup == (-1)) {
                      islandSize += (curBlob.size)
                    
                     if(edgeLabCount.toList.sum ==0){ //TODO Somehow this never occures. Maybe if we had an island on point 0,0,0 it would happen
                       remainingIslandBlobs.+=(altReturn._2)
                      // println(" [v:"+debug_countEdge+"]")
                  //     debug_countEdge=0
                     }
                     else{
                      
                      
                      val bestLabel = argmax(edgeLabCount) //TODO THIS IS WRONG. EDGELABCOUND can have 
                      //TODO if i want to include new clusters I need to also search for them here because centers is static
                      curBlob.foreach {
                        case (x, y, z) => {
                          blobConnections(x)(y)(z) = bestLabel
                          clusterAssign(x)(y)(z) = bestLabel
                        }
                      }
                      
                     }
                     
                 //    curBlob = new HashSet[(Int, Int, Int)]
                 //     edgeLabCount = Array.fill(centers.size) { 0 }

                    

                  } else if (foundGroup == clusterAssign(vX)(vY)(vZ)) {
                    curBlob.foreach(a => {
                      val x = a._1
                      val y = a._2
                      val z = a._3
                      blobConnections(x)(y)(z) = foundGroup  //TODO this is unnecessary alreayd done inside findConnection_f
                    })
                 //   curBlob = new HashSet[(Int, Int, Int)]
                  } else {
                    println("Error, should not get here #DASDWQ")
                  }
                }
              }
            }
          }
        }
      }

      //Connect Islands which did not find any neighbours (meaning they were somehow stuck inside another 
      
      def findFirstSetNeighbour (in: HashSet[(Int,Int,Int)]):Int={
        
        
        in.foreach( a =>{
          val x = a._1
          val y = a._2
          val z = a._3
          if(blobConnections(x+1)(y)(z)>=0)//TODO check bounds on all these 
            return blobConnections(x+1)(y)(z)
            if(blobConnections(x-1)(y)(z)>=0)
            return blobConnections(x-1)(y)(z)
            if(blobConnections(x)(y+1)(z)>=0)
            return blobConnections(x)(y+1)(z)
            if(blobConnections(x)(y-1)(z)>=0)
            return blobConnections(x)(y-1)(z)
            if(blobConnections(x)(y)(z+1)>=0)
            return blobConnections(x)(y)(z+1)
            if(blobConnections(x)(y)(z-1)>=0)
            return blobConnections(x)(y)(z-1)
            
            if(blobConnections(x+1)(y+1)(z)>=0)
            return blobConnections(x+1)(y+1)(z)
            if(blobConnections(x+1)(y-1)(z)>=0)
            return blobConnections(x+1)(y-1)(z)
            
            if(blobConnections(x+1)(y+1)(z+1)>=0)
            return blobConnections(x+1)(y+1)(z+1)
            if(blobConnections(x+1)(y-1)(z+1)>=0)
            return blobConnections(x+1)(y-1)(z+1)
            if(blobConnections(x+1)(y+1)(z-1)>=0)
            return blobConnections(x+1)(y+1)(z-1)
            if(blobConnections(x+1)(y-1)(z-1)>=0)
            return blobConnections(x+1)(y-1)(z-1)
            
            if(blobConnections(x+1)(y)(z+1)>=0)
            return blobConnections(x+1)(y)(z+1)
            if(blobConnections(x+1)(y)(z+1)>=0)
            return blobConnections(x+1)(y)(z+1)
            if(blobConnections(x+1)(y)(z-1)>=0)
            return blobConnections(x+1)(y)(z-1)
            if(blobConnections(x+1)(y)(z-1)>=0)
            return blobConnections(x+1)(y)(z-1)

        })
        
        
        return -1
        
      }
      
      /*
      val remainingOuter : scala.collection.parallel.mutable.ParSet[HashSet[(Int, Int, Int)]]= scala.collection.parallel.mutable.ParSet()
      remainingIslandBlobs.foreach( cordSet => {
      
      if(cordSet.size==0)
        println(" whySoZero")
        //Heuristic for speed is that we move in all 6 straight line directions of the average voxel to find the closest blob. 
        // But if we are an island inside an island 
        
        val foundGroup = findFirstSetNeighbour(cordSet)
        if(foundGroup>=0){
          cordSet.foreach(a => {
                      val x = a._1
                      val y = a._2
                      val z = a._3
                      blobConnections(x)(y)(z) = foundGroup  //TODO this is unnecessary alreayd done inside findConnection_f
                    })
                    
        }
        else{
          remainingOuter+=cordSet
        }
        
          
                      
      })
      print(" remainingOuter "+ remainingOuter.size+" ")
      */
      
      
      println("#Num Islands " + islandSize.size + " Avg Size(" + (islandSize.toList.sum / islandSize.size) + ")" + " max Size:(" + islandSize.toList.max + ")" + " min Size(" + islandSize.toList.min + ")"+" lostIslands "+remainingIslandBlobs.size)
      
    
    }

    var lastChange = Double.MaxValue
    var otherHalt = 0
    println(centers.size + " centers.size")
    while (otherHalt < 20) {

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
            if (!(vX < 0 || vX >= xDim || vY < 0 || vY >= yDim || vZ < 0 || vZ >= zDim)) {

              val curVox = Grey3dCord(vX, vY, vZ, aStack.getVoxel(vX, vY, vZ))
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

      val tCon = System.currentTimeMillis()
      enforceConnectedness()
      println((System.currentTimeMillis() - tCon) + " Enforce Connectedness")

      val updateCalcSum = Array.fill(centers.size) { Grey3dCord(0, 0, 0, 0.0) }
      val updateCalcCount = Array.fill(centers.size) { 0 }

      val tMore = System.currentTimeMillis()

      println("")
      println("How Many Unlabled Voxels per slice")
      (0 until xDim).toList.par.map { vX =>
        {
          var howManyNegOne = 0
          (0 until yDim).toList.map { vY =>
            {
              (0 until zDim).toList.map { vZ =>
                {

                  val usedForC = clusterAssign(vX)(vY)(vZ)
                  if (usedForC != (-1)) {
                    val myVoxVal = aStack.getVoxel(vX, vY, vZ)
                    val lastSum = updateCalcSum(usedForC)

                    updateCalcSum(usedForC) = Grey3dCord(vX + lastSum.x, vY + lastSum.y, vZ + lastSum.z, myVoxVal + lastSum.grey)
                    updateCalcCount(usedForC) = updateCalcCount(usedForC) + 1

                  } else {
                    howManyNegOne += 1
                  }
                }

              }
            }
          }
          print(" (" + howManyNegOne + ")")

        }
      }

      println((System.currentTimeMillis() - tMore) + " Update Cluster Centers")
      val lastCenters = centers.toArray
      var totalChange = 0.0
      val tTime = System.currentTimeMillis()
      var howManyCountZero = 0;
      for (cIDX <- 0 until centers.size) {
        val mySum = updateCalcSum(cIDX)
        val myCount = updateCalcCount(cIDX)
        if (myCount != 0) {
          centers(cIDX) = Grey3dCord(round(mySum.x / myCount),
            round(mySum.y / myCount),
            round(mySum.z / myCount),
            mySum.grey / myCount)
          totalChange += clusterDist(lastCenters(cIDX), cIDX)
        } else {
          howManyCountZero += 1
        }
      }

      println((System.currentTimeMillis() - tTime) + " Error Detection")

      println(" Centers Which where not used " + howManyCountZero)
      totalChange = totalChange / centers.size
      lastChange = totalChange
      printSuperPixels(label = "suppix_r_" + otherHalt)
      otherHalt += 1
      println(totalChange)

    }

    def printSuperPixels(borderCol: Double = 300.0, label: String = "NoName") {
      val imp2_pix = new Duplicator().run(imp2);
      val aStack_supPix = imp2_pix.getStack

      for (vX <- 0 until xDim) {
        for (vY <- 0 until yDim) {
          var lastLabel = -1
          for (vZ <- 0 until zDim) {
            if (clusterAssign(vX)(vY)(vZ) != lastLabel) {
              aStack_supPix.setVoxel(vX, vY, vZ, borderCol)
            }
            lastLabel = clusterAssign(vX)(vY)(vZ)
          }
        }
      }

      for (vX <- 0 until xDim) {
        for (vZ <- 0 until zDim) {
          var lastLabel = -1
          for (vY <- 0 until yDim) {

            if (clusterAssign(vX)(vY)(vZ) != lastLabel) {
              aStack_supPix.setVoxel(vX, vY, vZ, borderCol)
            }
            lastLabel = clusterAssign(vX)(vY)(vZ)
          }
        }
      }

      for (vZ <- 0 until zDim) {

        for (vY <- 0 until yDim) {
          var lastLabel = -1
          for (vX <- 0 until xDim) {
            if (clusterAssign(vX)(vY)(vZ) != lastLabel) {
              aStack_supPix.setVoxel(vX, vY, vZ, borderCol)
            }
            lastLabel = clusterAssign(vX)(vY)(vZ)
          }
        }
      }
      imp2_pix.setStack(aStack_supPix)
      //imp2_pix.show()
      IJ.saveAs(imp2_pix, "tif", "/home/mort/workspace/ScalaSLIC/data/image" + label + ".tif");
      //TODO free this memory after saving 
    }

    // ?? What is the gradiant in the 3x3 neighbourhood. For initialization of cluster centers (this is in he paper somewher)
    //Create a struct for saving the cluster centers  <color_stuff, x ,y ,z> 
    //Create a "hashFunc" which maps pixels to a spot in the grid. 
    //Write distance measure which can be easily adapted based on the type of image data
    //Pixel to cluster assignment matrix  same size as original 
    //MinDistance matrix for each voxel store the distance to its nearest cluster
    // ?? What is the residual error in this case. Assign pixels to their cluster and see the diff i guess. 

  }

}