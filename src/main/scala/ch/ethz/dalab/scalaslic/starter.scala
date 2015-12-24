package ch.ethz.dalab.scalaslic

import scala.util.Try
import ij.io.Opener
import java.io.File

object starter {

 

  def main(args: Array[String]): Unit = {

    if (args.length < 2) {
      println("Incorrect number of arguments"); 
      throw new Exception;
    }

    val fileLoc = args(0);
    val superpixelSize = args(1).toInt;  
    val compactness = ((args.lift(2)).getOrElse(Double.MinValue.toString())).toDouble
    val minSuperPixelSize = (args.lift(3)).getOrElse("-1").toInt
    val imageData:Array[Array[Array[Int]]] =utils.fileReadling(fileLoc);
    
    val slicworker = new SLICgreySimple(imageData,superpixelSize,compactness,myMinSize=minSuperPixelSize);
    val superPixelMask= slicworker.calcSuperPixels();
    
    
    val opener = new Opener();
    val origFile = opener.openImage(fileLoc);
    println(new File(".").getAbsolutePath()); 
    
    utils.printSuperPixels(superPixelMask,origFile,300.0,new File(fileLoc).getName()+"(S_"+superpixelSize+") (M_"+compactness+") (MinSupSiz_"+minSuperPixelSize+")");
    println("Done.");
  }
}