use std::vec::Vec;
use std::fs::File;
use std::io::{BufWriter, Write, Error};
use color::Color;

pub struct Image {
  width: usize,
  height: usize,
  pixels: Vec<Color>,
}

fn to_integer<T: Into<f64>>(x: T) -> u64 {
  let mut x = x.into();

  if x < 0.0 { x = 0.0; } else
  if x > 1.0 { x = 1.0; };

  // Apply gamma correction and convert to integer.
  (x.powf(1.0 / 2.2) * 255.0 + 0.5) as u64
}

impl Image {
  pub fn new(width: usize, height: usize) -> Image {
    let size = (width * height) as usize;
    Image { width: width, height: height, pixels: vec![Color::new(0.0, 0.0, 0.0); size] }
  }

  pub fn get_pixel_mut(&mut self, x: usize, y: usize) -> &mut Color {
    let image_index = (self.height - y - 1) * self.width + x;
    &mut self.pixels[image_index]
  }

  pub fn save(&self, file_name: &str) -> Result<File, Error> {
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
