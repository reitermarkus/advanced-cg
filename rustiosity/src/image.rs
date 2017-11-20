use std::vec::Vec;
use std::fs::File;
use std::io::{BufWriter, Write, Error};
use color::Color;

pub struct Image {
  width: i64,
  height: i64,
  pixels: Vec<Color>,
}

fn to_integer<T: Into<f64>>(x: T) -> usize {
  let mut x = x.into();

  if x < 0.0 { x = 0.0; } else
  if x > 1.0 { x = 1.0; };

  // Apply gamma correction and convert to integer.
  (x.powf(1.0 / 2.2) * 255.0 + 0.5) as usize
}

impl Image {
  pub fn new(width: i64, height: i64) -> Image {
    let size = (width * height) as usize;
    Image { width: width, height: height, pixels: vec![Color::new(0.0, 0.0, 0.0); size] }
  }

  fn index(&self, x: i64, y: i64) -> usize {
    ((self.height - y - 1) * self.width + x) as usize
  }

  pub fn add_color(&mut self, x: i64, y: i64, color: &Color) {
    let image_index = self.index(x, y);
    self.pixels[image_index] = self.pixels[image_index] + color;
  }

  pub fn save(&self, file_name: &String) -> Result<File, Error> {
    let file = File::create(file_name)?;

    {
      let mut writer = BufWriter::new(&file);

      write!(writer, "P3\n{} {}\n{}\n", self.width, self.height, 255)?;

      for pixel in &self.pixels {
        let r = to_integer(pixel.x);
        let g = to_integer(pixel.y);
        let b = to_integer(pixel.z);
        write!(writer, "{} {} {} ", r, g, b)?;
      }
    }

    Ok(file)
  }
}
