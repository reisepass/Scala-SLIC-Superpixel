name := "ScalaSLIC"

organization := "ch.ethz.dalab"

version := "0.1-SNAPSHOT"

scalaVersion := "2.10.4"
javaOptions in run ++= Seq("-Xms256M", "-Xmx2G", "-XX:MaxPermSize=1024M", "-XX:+UseConcMarkSweepGC")

libraryDependencies += "gov.nih.imagej" % "imagej" % "1.46"
libraryDependencies += "org.scalanlp" %% "breeze" % "0.11.2"


resolvers += Resolver.sonatypeRepo("public")

