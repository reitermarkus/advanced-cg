mod color;
mod vector;
mod triangular;
mod triangle;
mod image;
mod patch_triangle;
mod ray;

extern crate crossbeam;
extern crate ordered_float;
extern crate rand;
extern crate rayon;

use rayon::prelude::*;
use rand::distributions::{IndependentSample, Range};

use std::vec::Vec;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::io::{stdout, Write};
use std::time::{Duration, Instant};
use std::sync::{Arc};
use std::thread;

use image::Image;
use ray::Ray;
use triangular::Triangular;
use triangle::Triangle;
use patch_triangle::PatchTriangle;
use vector::Vector;
use color::Color;

fn tris() -> Vec<Triangle> {
  vec![
  /* Cornell Box walls */
    Triangle::new(Vector::new(  0.0,  0.0,   0.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Back:   bottom-left
    Triangle::new(Vector::new(100.0, 80.0,   0.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Back:   top-right
    Triangle::new(Vector::new(  0.0,  0.0, 170.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,   0.0, -170.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Bottom: front-left
    Triangle::new(Vector::new(100.0,  0.0,   0.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,   0.0,  170.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Bottom: back-right
    Triangle::new(Vector::new(  0.0, 80.0,   0.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,   0.0,  170.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Top:    back-left
    Triangle::new(Vector::new(100.0, 80.0, 170.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,   0.0, -170.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Top:    front-right
    Triangle::new(Vector::new(  0.0,  0.0, 170.0), Vector::new(   0.0, 0.0, -170.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.75, 0.25, 0.25)), // Left:   front-bottom
    Triangle::new(Vector::new(  0.0, 80.0,   0.0), Vector::new(   0.0, 0.0,  170.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.75, 0.25, 0.25)), // Left:   back-top
    Triangle::new(Vector::new(100.0,  0.0,   0.0), Vector::new(   0.0, 0.0,  170.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.25, 0.25, 0.75)), // Right:  back-bottom
    Triangle::new(Vector::new(100.0, 80.0, 170.0), Vector::new(   0.0, 0.0, -170.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.25, 0.25, 0.75)), // Right:  front-top
    Triangle::new(Vector::new(100.0,  0.0, 170.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.25, 0.75,  0.25)),  // Front:  bottom-right (not visible)
    Triangle::new(Vector::new(  0.0, 80.0, 170.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.25, 0.75,  0.25)),  // Front:  top-left (not visible)

    /* Area light source on top */
    Triangle::new(Vector::new(40.0, 79.99, 65.0), Vector::new( 20.0, 0.0, 0.0), Vector::new(0.0, 0.0,  20.0), Color::new(12.0, 12.0, 12.0), Color::new(0.75, 0.75, 0.75)), // back-left
    Triangle::new(Vector::new(60.0, 79.99, 85.0), Vector::new(-20.0, 0.0, 0.0), Vector::new(0.0, 0.0, -20.0), Color::new(12.0, 12.0, 12.0), Color::new(0.75, 0.75, 0.75)), // front-right

    /* Cuboid in room */
    Triangle::new(Vector::new(30.0,  0.0, 100.0), Vector::new(  0.0, 0.0, -20.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Right: front-bottom
    Triangle::new(Vector::new(30.0, 40.0,  80.0), Vector::new(  0.0, 0.0,  20.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Right: back-top
    Triangle::new(Vector::new(10.0,  0.0,  80.0), Vector::new(  0.0, 0.0,  20.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Left:  back-bottom
    Triangle::new(Vector::new(10.0, 40.0, 100.0), Vector::new(  0.0, 0.0, -20.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Left:  front-top
    Triangle::new(Vector::new(10.0,  0.0, 100.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Front: bottom-left
    Triangle::new(Vector::new(30.0, 40.0, 100.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Front: top-right
    Triangle::new(Vector::new(30.0,  0.0,  80.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Back:  bottom-right
    Triangle::new(Vector::new(10.0, 40.0,  80.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Back:  top-left
    Triangle::new(Vector::new(10.0, 40.0, 100.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0,   0.0, -20.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Top:   front-left
    Triangle::new(Vector::new(30.0, 40.0,  80.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0,   0.0,  20.0), Color::zero(), Color::new(0.75, 0.75, 0.75)), // Top:   back-right
  ]
}

fn intersect_scene(tris: &[Triangle], ray: &Ray, t: &mut f64, id: &mut i64, normal: &mut Vector) -> bool {
  *t = 1e20;
  *id = -1;

  for i in 0..tris.len() {
    let d = tris[i].intersect(ray);
    if d > 0.0 && d < *t {
      *t = d;
      *id = i as i64;
      *normal = tris[i].normal;
    }
  }

  *t < 1e20
}

fn calculate_form_factors(tris: &mut [Triangle], divisions: u64, mc_sample: i64) -> HashMap<usize, Vec<HashMap<usize, Vec<f64>>>> {
  let mut patch_num = 0;

  // Total number of patches in scene.
  let n = tris.len();
  for i in 0..n {
    tris[i].init_patches(divisions);
    patch_num += tris[i].patches.len();
  }

  let form_factor_num = patch_num.pow(2);

  println!("Number of triangles: {}", n);
  println!("Number of patches: {}", patch_num);
  println!("Number of form factors: {}", form_factor_num);

  // Loop over all triangles in scene.
  let mut form_factors: HashMap<_, _> = (0..n).into_par_iter().map(|i| {
    let mut maps: Vec<HashMap<usize, Vec<f64>>> = vec![HashMap::new(); tris[i].patches.len()];

    // Loop over all patches in rectangle i.
    for p_i in 0..tris[i].patches.len() {
      for j in 0..n {
        maps[p_i].insert(j, vec![0.0; tris[j].patches.len()]);
      }

      // Loop over all triangles in scene for triangles i.
      for j in (i + 1)..n {
        // Loop over all patches in rectangle j.
        for p_j in 0..tris[j].patches.len() {
          // Monte Carlo integration of form factor double integral.

          let t_i: &PatchTriangle = &tris[i].patches[p_i];
          let t_j: &PatchTriangle = &tris[j].patches[p_j];

          // Uniform PDF for Monte Carlo (1 / Ai) x (1 / Aj).
          let pdf = 1.0 / (t_i.area * t_j.area);

          // Determine rays of NixNi uniform samples of patch
          // on i to NjxNj uniform samples of patch on j.
          let mut form_factor = (0..mc_sample).into_iter().map(|_| {
            let xi: Vector = t_i.random_sample();
            let xj: Vector = t_j.random_sample();

            // Check for visibility between sample points.
            let ij: Vector = (&xj - &xi).normalize();

            let mut t = 0.0;
            let mut id = -1;
            let mut normal = Vector::zero();
            if intersect_scene(&tris, &Ray::new(&xi, &ij), &mut t, &mut id, &mut normal) && id != j as i64 {
              return 0.0; // If intersection with other triangle.
            }

            // Cosines of angles beteen normals and ray inbetween.
            let d0 = tris[i].normal.dot_product(&ij);
            let d1 = tris[j].normal.dot_product(&(-1.0 * ij));

            // Continue if patches facing each other.
            if d0 <= 0.0 || d1 <= 0.0 {
              return 0.0;
            }

            // Sample form factor.
            d0 * d1 / (PI * (&xj - &xi).length_squared())
          }).sum();

          form_factor /= pdf;

          // Divide by number of samples.
          form_factor /= mc_sample as f64;

          maps[p_i].get_mut(&j).unwrap().insert(p_j, form_factor);
        }
      }
    }

    print!("*");
    stdout().flush().unwrap();

    (i, maps)
  }).collect();

  println!();

  for i in 0..tris.len() {
    for p_i in 0..tris[i].patches.len() {
      for j in (i + 1)..tris.len() {
        for p_j in 0..tris[j].patches.len() {
          form_factors.get_mut(&j).unwrap()[p_j].get_mut(&i).unwrap()[p_i] = form_factors.get(&i).unwrap()[p_i].get(&j).unwrap()[p_j];
        }
      }
    }
  }

  // Divide by area to get final form factors.
  for i in 0..tris.len() {
    for p_i in 0..tris[i].patches.len() {
      for j in 0..tris.len() {
        for p_j in 0..tris[j].patches.len() {
          let area = tris[i].patches[p_i].area;
          let form_factor = form_factors[&i][p_i][&j][p_j] / area;
          form_factors.get_mut(&i).unwrap()[p_i].get_mut(&j).unwrap()[p_j] = if form_factor < 0.0 { 0.0 } else if form_factor > 1.0 { 1.0 } else { form_factor};
        }
      }
    }
  }

  form_factors
}

fn calculate_radiosity(tris: &mut [Triangle], form_factors: &HashMap<usize, Vec<HashMap<usize, Vec<f64>>>>) {
  (0..tris.len()).into_par_iter().for_each(|i| {
    for p_a in 0..tris[i].patches.len() {
      let color = (0..tris.len()).map(|j| {
        if i == j { return Color::zero(); }
        (0..tris[j].patches.len()).map(|p_b| {
          form_factors[&i][p_a][&j][p_b] * tris[j].patches[p_b].get_color()
        }).sum()
      }).sum();

      // Multiply sum with color of patch and add emission and
      // store overall patch radiosity of current iteration.
      tris[i].patches[p_a].set_color(tris[i].color.entrywise_product(&color) + tris[i].emission);
    }
  });
}

fn calculate_vertex_colors(tris: &[Triangle]) -> HashMap<Vector, HashMap<Vector, Color>> {
  let mut vertex_colors: HashMap<Vector, HashMap<Vector, Color>> = HashMap::new();
  let mut vertex_counts: HashMap<Vector, HashMap<Vector, u64>> = HashMap::new();

  for tri in tris {
    vertex_colors.insert(tri.normal, HashMap::with_capacity(tri.patches.len()));
    vertex_counts.insert(tri.normal, HashMap::with_capacity(tri.patches.len()));
  }

  for tri in tris {
    for p in 0..tri.patches.len() {
      let vertex_counts = vertex_counts.get_mut(&tri.normal).unwrap();
      let vertex_colors = vertex_colors.get_mut(&tri.normal).unwrap();

      vertex_colors.entry(tri.patches[p].a).or_insert(Color::zero());
      vertex_colors.entry(tri.patches[p].b).or_insert(Color::zero());
      vertex_colors.entry(tri.patches[p].c).or_insert(Color::zero());
      vertex_counts.entry(tri.patches[p].a).or_insert(0);
      vertex_counts.entry(tri.patches[p].b).or_insert(0);
      vertex_counts.entry(tri.patches[p].c).or_insert(0);

      *vertex_counts.get_mut(&tri.patches[p].a).unwrap() += 1;
      *vertex_counts.get_mut(&tri.patches[p].b).unwrap() += 1;
      *vertex_counts.get_mut(&tri.patches[p].c).unwrap() += 1;
      *vertex_colors.get_mut(&tri.patches[p].a).unwrap() += tri.patches[p].get_color();
      *vertex_colors.get_mut(&tri.patches[p].b).unwrap() += tri.patches[p].get_color();
      *vertex_colors.get_mut(&tri.patches[p].c).unwrap() += tri.patches[p].get_color();
    }
  }

  for (normal, colors) in vertex_colors.iter_mut() {
    for (vertex, color) in colors.iter_mut() {
      *color /= *vertex_counts.get(&normal).unwrap().get(&vertex).unwrap() as f64;
    }
  }

  return vertex_colors;
}

fn radiance(tris: &[Triangle], ray: &Ray, vertex_colors: &HashMap<Vector, HashMap<Vector, Color>>) -> (Color, Color) {
  let background_color = Color::zero();

  let mut t = -1.0;
  let mut id = -1;
  let mut normal = Vector::zero();

  // Find intersected triangle.
  if !intersect_scene(tris, ray, &mut t, &mut id, &mut normal) {
    return (background_color, background_color);
  }

  // Find intersected patch.
  let obj: &Triangle = &tris[id as usize];

  let find_patch_index = |triangle: &Triangle, ray: &Ray| -> Option<usize> {
    for p in 0..triangle.patches.len() {
      if triangle.patches[p].intersect(ray) != 0.0 {
        return Some(p);
      }
    }

    None
  };

  let idx = match find_patch_index(obj, ray) {
    Some(p) => p,
    None => panic!("Could not find index of patch."),
  };

  let hitpoint = ray.org + t * ray.dir;

  let bary = obj.patches[idx].barycentric_coordinates_at(&hitpoint);

  let a: Color = vertex_colors[&obj.normal][&obj.patches[idx].a];
  let b: Color = vertex_colors[&obj.normal][&obj.patches[idx].b];
  let c: Color = vertex_colors[&obj.normal][&obj.patches[idx].c];

  let interpolated: Color = a * bary.x.into_inner() + b * bary.z.into_inner() + c * bary.y.into_inner();

  return (obj.patches[idx].get_color(), interpolated);
}


fn main() {
  let total_bench = Instant::now();
  let mut tris = tris();

  let divisions: u64 = 10;
  let mc_samples = 10;

  println!("Calculating form factors ...");
  let form_factor_bench = Instant::now();
  let form_factors = calculate_form_factors(&mut tris, divisions, mc_samples);
  let form_factor_bench_elapsed = form_factor_bench.elapsed();

  println!("Calculating radiosity ...");
  let iterations = 40;
  let radiosity_bench = Instant::now();
  for i in 0..iterations {
    print!("{} ", i);
    calculate_radiosity(&mut tris, &form_factors);
  }
  let radiosity_bench_elapsed = radiosity_bench.elapsed();
  println!();

  println!("Calculating vertex colors ...");
  let vertex_colors = calculate_vertex_colors(&tris);

  println!("Rendering and saving images ...");

  let width = 640;
  let height = 480;

  let samples = 4;

  let camera = Ray::new(&Vector::new(50.0, 52.0, 295.6), &Vector::new(0.0, -0.042612, -1.0).normalize());

  // Image edge vectors for pixel sampling.
  let cx: Vector = Vector::new(width as f64 * 0.5135 / height as f64, 0.0, 0.0);
  let cy: Vector = cx.cross_product(&camera.dir).normalize() * 0.5135;

  let image = Arc::new(Image::new(width, height));
  let image_interpolated = Arc::new(Image::new(width, height));

  let image_clone = image.clone();
  let image_interpolated_clone = image_interpolated.clone();

  let rendering_and_saving_bench = Instant::now();

  let image_thread = thread::spawn(move || {
    thread::park();

    if let Err(e) = image_clone.save("image_patches.ppm") {
      panic!("{:?}", e)
    }
  });

  let image_interpolated_thread = thread::spawn(move || {
    thread::park();

    if let Err(e) = image_interpolated_clone.save("image_smooth.ppm") {
      panic!("{:?}", e)
    }
  });

  let start_saving = if height < 20 { height / 2 } else { height - 20 };

  // Loop over image rows.
  (0..height).into_par_iter().for_each(|y| {
    let between = Range::new(0.0, 1.0);
    let mut rng = rand::thread_rng();

    // Loop over image columns.
    for x in 0..width {
      let mut accumulated_radiance = Color::zero();
      let mut accumulated_radiance_interpolated = Color::zero();

      // 2 x 2 subsampling per pixel.
      for sy in 0..2 {
        for sx in 0..2 {

          // Computes radiance at subpixel using multiple samples.
          for _ in 0..samples {
            let mut nu_filter_samples = || -> f64 {
              let r = 2.0 * between.ind_sample(&mut rng) as f64;
              if r < 1.0 { r.sqrt() - 1.0 } else { 1.0 - (2.0 - r).sqrt() }
            };

            let dx = nu_filter_samples();
            let dy = nu_filter_samples();

            // Ray direction into scene from camera through sample.
            let dir: Vector = cx * (((x as f64) + ((sx as f64) + 0.5 + dx) / 2.0) / (width as f64) - 0.5) +
                              cy * ((((height - y - 1) as f64) + ((sy as f64) + 0.5 + dy) / 2.0) / (height as f64) - 0.5) +
                              camera.dir;

            // Extend camera ray to start inside box.
            let start: Vector = camera.org + dir * 130.0;

            let (color, color_interpolated) = radiance(&tris, &Ray::new(&start, &dir.normalize()), &vertex_colors);

            // Determine constant radiance.
            accumulated_radiance += color;

            // Determine interpolated radiance.
            accumulated_radiance_interpolated += color_interpolated;
          }
        }
      }

      image.set_color(x, y, accumulated_radiance / samples as f64);
      image_interpolated.set_color(x, y, accumulated_radiance_interpolated / samples as f64);

      if y == start_saving {
        image_thread.thread().unpark();
        image_interpolated_thread.thread().unpark();
      }
    }
  });

  image_thread.join().unwrap();
  image_interpolated_thread.join().unwrap();

  let rendering_and_saving_bench_elapsed = rendering_and_saving_bench.elapsed();
  let total_bench_elapsed = total_bench.elapsed();

  let into_ms = |x: Duration| -> u64 {
    return (x.as_secs() * 1_000) + (x.subsec_nanos() / 1_000_000) as u64;
  };

  println!("┌────────────────────┬───────────┐");
  println!("│ Form Factors       │ {:#6 } ms │", into_ms(form_factor_bench_elapsed));
  println!("├────────────────────┼───────────┤");
  println!("│ Radiosity          │ {:#6 } ms │", into_ms(radiosity_bench_elapsed));
  println!("├────────────────────┼───────────┤");
  println!("│ Rendering & Saving │ {:#6 } ms │", into_ms(rendering_and_saving_bench_elapsed));
  println!("┢━━━━━━━━━━━━━━━━━━━━╈━━━━━━━━━━━┪");
  println!("┃ Total              ┃ {:#6 } ms ┃", into_ms(total_bench_elapsed));
  println!("┗━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━┛");
}
