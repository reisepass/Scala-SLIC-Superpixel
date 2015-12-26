# Scala-SLIC-Superpixel
Simple Linear Iterative Clustering Superpixel Algorithm in Scala

Example Usage for greyscale images. 

$ sbt "run-main ch.ethz.dalab.scalaslic.starter exampleData/LennaGray.png 30 60 -1"


The system is easily extendable to color or any othe grid datastructure by simply defining a distance function and a rolling Average function.



The Original SLIC algorithm : 
  http://ivrl.epfl.ch/research/superpixels 


First published: 
  Radhakrishna Achanta, Appu Shaji, Kevin Smith, 
  Aurelien Lucchi, Pascal Fua, and Sabine Süsstrunk,
  SLIC Superpixels, EPFL Technical Report no. 149300, June 2010.
