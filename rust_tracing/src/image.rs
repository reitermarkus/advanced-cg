use std::fs::File;
use std::io::{BufWriter, Write, Error};
use std::sync::atomic::Ordering::{Acquire, Release};

use crossbeam::epoch::{self, Atomic, Owned, CompareAndSetError};

use color::Color;

const GAMMA: f64 = 2.2;

pub struct Image {
  width: usize,
  height: usize,
  pixels: Vec<Atomic<Color>>,
}

#[inline(always)]
fn to_byte<T: Into<f64>>(x: T) -> u8 {
  let mut x = x.into();

  if x < 0.0 { x = 0.0; } else
  if x > 1.0 { x = 1.0; };

  // Apply gamma correction and convert to integer.
  (x.powf(1.0 / GAMMA) * 255.0 + 0.5) as u8
}

impl Image {
  pub fn new(width: usize, height: usize) -> Image {
    let size = (width * height) as usize;

    Image { width: width, height: height, pixels: (0..size).map(|_| Atomic::null()).collect() }
  }

  fn index(&self, x: usize, y: usize) -> usize {
    y * self.width + x
  }

  pub fn set_color(&self, x: usize, y: usize, color: Color) {
    let image_index = self.index(x, y);
    let mut color = Owned::new(color);

    let guard = epoch::pin();

    loop {
      let c = self.pixels[image_index].load(Acquire, &guard);

      match self.pixels[image_index].compare_and_set(c, color, Release, &guard) {
        Ok(_) => return,
        Err(CompareAndSetError { new, .. }) => color = new,
      }
    }
  }

  pub fn save(&self, file_name: &str) -> Result<File, Error> {
    let file = File::create(file_name)?;

    {
      let mut writer = BufWriter::new(&file);

      write!(writer, "P6\n{} {}\n{}\n", self.width, self.height, 255)?;

      for i in 0..self.pixels.len() {
        let guard = epoch::pin();

        loop {
          let color = self.pixels[i].load(Acquire, &guard);

          if color.is_null() {
            continue
          }

          let Color { x: r, y: g, z: b } = *unsafe { color.deref() };

          writer.write(&[to_byte(r), to_byte(g), to_byte(b)])?;

          break;
        }
      }
    }

    Ok(file)
  }
}

unsafe impl Sync for Image {}
unsafe impl Send for Image {}
