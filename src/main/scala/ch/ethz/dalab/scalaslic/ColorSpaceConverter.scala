package ch.ethz.dalab.scalaslic

object ColorSpaceConverter {
  case class RGB(r: Double, g: Double, b: Double)
  case class LAB(l: Double, a: Double, b: Double)
  case class HSV(h: Double, s: Double, v: Double)

  def rgbToLab(rgb: RGB): LAB = {
    // First convert RGB to XYZ
    val r = rgb.r / 255.0
    val g = rgb.g / 255.0
    val b = rgb.b / 255.0

    // Convert to XYZ
    val x = r * 0.4124564 + g * 0.3575761 + b * 0.1804375
    val y = r * 0.2126729 + g * 0.7151522 + b * 0.0721750
    val z = r * 0.0193339 + g * 0.1191920 + b * 0.9503041

    // Convert XYZ to LAB
    val xn = 95.047
    val yn = 100.000
    val zn = 108.883

    val fx = if (x / xn > 0.008856) Math.pow(x / xn, 1.0/3.0) else (7.787 * x / xn + 16.0/116.0)
    val fy = if (y / yn > 0.008856) Math.pow(y / yn, 1.0/3.0) else (7.787 * y / yn + 16.0/116.0)
    val fz = if (z / zn > 0.008856) Math.pow(z / zn, 1.0/3.0) else (7.787 * z / zn + 16.0/116.0)

    val l = 116.0 * fy - 16.0
    val a = 500.0 * (fx - fy)
    val b = 200.0 * (fy - fz)

    LAB(l, a, b)
  }

  def rgbToHsv(rgb: RGB): HSV = {
    val r = rgb.r / 255.0
    val g = rgb.g / 255.0
    val b = rgb.b / 255.0
    
    val max = math.max(math.max(r, g), b)
    val min = math.min(math.min(r, g), b)
    val delta = max - min
    
    val h = if (delta == 0) 0
            else if (max == r) 60 * ((g - b) / delta % 6)
            else if (max == g) 60 * ((b - r) / delta + 2)
            else 60 * ((r - g) / delta + 4)
    
    val s = if (max == 0) 0 else delta / max
    val v = max
    
    HSV(h, s, v)
  }
} 