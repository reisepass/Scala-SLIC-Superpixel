
import breeze.linalg._
import breeze.stats.DescriptiveStats._
import breeze.stats._
import breeze.numerics._
import scala.collection.mutable.HashSet
import scala.collection.mutable.ListBuffer
import scala.util.Random
import scala.collection.parallel.mutable.ParHashSet
import scala.collection.parallel.mutable.ParMap


/*
 *  DataType for simple greyscale images will just be Double. for RGB it could be Int or (Int,Int,Int)  for example 
 *  dataDistanceFn just needs to define the distance between two colors. For grey scale we use linear euclidian distnace
 *  S is the grid size for the initial placement of super pixel centers. 
 *  
 */
class SLIC[DataType](distFn:(DataType,DataType)=>Double, 
    sumFn:((DataType,DataType)=>DataType), 
    normFn:((DataType,Int)=>DataType),
    image: Array[Array[Array[DataType]]], S:Int , debug:Boolean=true) {


//TODO try out different tif image stacks and see what dim2 is 
    
    val xDim = image.size
    val yDim = image(0).size
    val zDim = image(0)(0).size 
    if(debug)
    println("(" + xDim + "," + yDim + "," + zDim + ") image Dims")
    assert(xDim > 0 & yDim > 0 & zDim > 0)
    case class DatumCord[DataType](x:Int,y:Int,z:Int,cont:DataType)
    
    def placeInitSuperCenters(gridInterval: Int): Array[(Int, Int, Int)] = {
      val out = for {
        x_s <- round(gridInterval / 2) until xDim by gridInterval;
        y_s <- round(gridInterval / 2) until yDim by gridInterval;
        z_s <- round(gridInterval / 2) until zDim by gridInterval
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
            
            difSum += distFn(myCol,image(dx+x+1)(dy+y)(dz+z))
            difSum += distFn(myCol,image(dx+x-1)(dy+y)(dz+z))
            difSum += distFn(myCol,image(dx+x)(dy+y+1)(dz+z))
            difSum += distFn (myCol,image(dx+x)(dy+y-1)(dz+z))
            difSum += distFn (myCol,image(dx+x)(dy+y)(dz+z+1))
            difSum += distFn (myCol,image(dx+x)(dy+y)(dz+z-1))
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
    val centerCoords = moveInitSuperCentersToLowContrast(initCenter, 3)
    val clusterAssign = Array.fill(xDim, yDim, zDim) { -1 }.par //TODO Can be changed to Range for less memory usage
    println(" clusterAssign size = " + clusterAssign.size + clusterAssign.toString())
    val lastDistance = Array.fill(xDim, yDim, zDim) { Double.MaxValue }.par
    

    val centers = centerCoords.map(a => DatumCord(a._1, a._2, a._3, image(a._1)(a._2)(a._3))).par 

    val t0 = System.currentTimeMillis()
    println((t0 - tt0) + " 58 to 66")
    val clusterMaxColDist = Array.fill(centers.size) { Double.MinValue }
    val t1 = System.currentTimeMillis()
    println((t0 - t1) + " for clusterMax Inst")
    val clusterMaxSpaceDist = Array.fill(centers.size) { Double.MinValue }
    println((System.currentTimeMillis() - t1) + "line 75")
    
     def clusterDist(point: DatumCord[DataType], centerID: Int): Double = {
      val center = centers(centerID)
      val d_c = distFn(center.cont,point.cont)
      val d_s = sqrt(Math.pow(point.x - center.x, 2) + Math.pow(point.y - center.y, 2) + Math.pow(point.z - center.z, 2))
      if (clusterMaxColDist(centerID) < d_c)
        clusterMaxColDist(centerID) = d_c
      if (clusterMaxSpaceDist(centerID) < d_s)
        clusterMaxSpaceDist(centerID) = d_s
      //TODO question: should i update the clusterMax all the time or just after each round ? 
      sqrt(Math.pow(d_c / clusterMaxColDist(centerID), 2) + Math.pow(d_s / clusterMaxSpaceDist(centerID), 2))

    }


    def enforceConnectedness() = {
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
      val newCenters = new ListBuffer[(DatumCord[DataType], Int)]

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

      val tCon = System.currentTimeMillis()
      enforceConnectedness()
      println((System.currentTimeMillis() - tCon) + " Enforce Connectedness")


      val updateCalcSum = ParMap[Int,DatumCord[DataType]]()
      val updateCalcCount = ParMap[Int,Int]()

      val tMore = System.currentTimeMillis()

      println("")
      println("How Many Unlabled Voxels per slice")
     (0 until xDim).toList.map { vX => //TODO parallize this.  (ParMap seems to not syncronize the way i expected it...)  
        {
          var howManyNegOne = 0
          (0 until yDim).toList.map { vY =>
            {
              (0 until zDim).toList.map { vZ =>
                {

                  val usedForC = clusterAssign(vX)(vY)(vZ)
                  if (usedForC != (-1)) {
                    if(updateCalcSum.contains(usedForC)){
                      val old  =updateCalcSum(usedForC)
                      updateCalcSum.put(usedForC, DatumCord(vX+old.x,vY+old.y,vZ+old.z,sumFn(image(vX)(vY)(vZ),old.cont)))
                      updateCalcCount.put(usedForC, updateCalcCount(usedForC)+1)
                    }
                    else{
                      updateCalcSum.put(usedForC, DatumCord(vX,vY,vZ,image(vX)(vY)(vZ)))
                      updateCalcCount.put(usedForC,1)
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
      
      updateCalcSum.keys.foreach { key =>{
        val old = updateCalcSum(key)
        val count = updateCalcCount(key)
         updateCalcSum.put(key,DatumCord(round(old.x/count),round(old.y/count),round(old.z/count),normFn(updateCalcSum(key).cont,count)))
         }
      }
      
      

      println((System.currentTimeMillis() - tMore) + " Update Cluster Centers")
      val lastCenters = centers.toArray
      var totalChange = 0.0
      val tTime = System.currentTimeMillis()
      var howManyCountZero = 0;
      for (cIDX <- 0 until centers.size) {
        if(updateCalcCount.contains(cIDX)){
          totalChange += clusterDist(lastCenters(cIDX), cIDX)
        } else {
          howManyCountZero += 1
        }
      }

      println((System.currentTimeMillis() - tTime) + " Error Detection")

      println(" Centers Which where not used " + howManyCountZero)
      totalChange = totalChange / centers.size
      lastChange = totalChange
    
      otherHalt += 1
      println(totalChange)

    
    
    
    }
 
    def clusterAssign(x:Int,y:Int,z:Int):Int=clusterAssign(x)(y)(z)
    
 

}