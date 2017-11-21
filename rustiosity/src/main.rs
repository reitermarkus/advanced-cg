mod color;
mod vector;
mod triangular;
mod triangle;
mod image;
mod patch_triangle;
mod ray;

extern crate rayon;
extern crate rand;

use rayon::prelude::*;
use rand::distributions::{IndependentSample, Range};

use std::vec::Vec;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::io::{self, Write};
use std::sync::{mpsc, Mutex};
use std::time::Instant;

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
    Triangle::new(Vector::new(  0.0,  0.0,   0.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,  80.0,    0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Back:   bottom-left
    Triangle::new(Vector::new(100.0, 80.0,   0.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0, -80.0,    0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Back:   top-right
    Triangle::new(Vector::new(  0.0,  0.0, 170.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,   0.0, -170.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Bottom: front-left
    Triangle::new(Vector::new(100.0,  0.0,   0.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,   0.0,  170.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Bottom: back-right
    Triangle::new(Vector::new(  0.0, 80.0,   0.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,   0.0,  170.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Top:    back-left
    Triangle::new(Vector::new(100.0, 80.0, 170.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,   0.0, -170.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Top:    front-right
    Triangle::new(Vector::new(  0.0,  0.0, 170.0), Vector::new(   0.0, 0.0, -170.0), Vector::new(0.0,  80.0,    0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.25, 0.25)), // Left:   front-bottom
    Triangle::new(Vector::new(  0.0, 80.0,   0.0), Vector::new(   0.0, 0.0,  170.0), Vector::new(0.0, -80.0,    0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.25, 0.25)), // Left:   back-top
    Triangle::new(Vector::new(100.0,  0.0,   0.0), Vector::new(   0.0, 0.0,  170.0), Vector::new(0.0,  80.0,    0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.25, 0.25, 0.75)), // Right:  back-bottom
    Triangle::new(Vector::new(100.0, 80.0, 170.0), Vector::new(   0.0, 0.0, -170.0), Vector::new(0.0, -80.0,    0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.25, 0.25, 0.75)), // Right:  front-top
    Triangle::new(Vector::new(100.0,  0.0, 170.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,  80.0,    0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.25, 0.75,  0.25)),  // Front:  bottom-right (not visible)
    Triangle::new(Vector::new(  0.0, 80.0, 170.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0, -80.0,    0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.25, 0.75,  0.25)),  // Front:  top-left (not visible)

    /* Area light source on top */
    Triangle::new(Vector::new(40.0, 79.99, 65.0), Vector::new( 20.0, 0.0, 0.0), Vector::new(0.0, 0.0,  20.0), Color::new(12.0, 12.0, 12.0), Color::new(0.75, 0.75, 0.75)), // back-left
    Triangle::new(Vector::new(60.0, 79.99, 85.0), Vector::new(-20.0, 0.0, 0.0), Vector::new(0.0, 0.0, -20.0), Color::new(12.0, 12.0, 12.0), Color::new(0.75, 0.75, 0.75)), // front-right

    /* Cuboid in room */
    Triangle::new(Vector::new(30.0,  0.0, 100.0), Vector::new(  0.0, 0.0, -20.0), Vector::new(0.0,  40.0,   0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Right: front-bottom
    Triangle::new(Vector::new(30.0, 40.0,  80.0), Vector::new(  0.0, 0.0,  20.0), Vector::new(0.0, -40.0,   0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Right: back-top
    Triangle::new(Vector::new(10.0,  0.0,  80.0), Vector::new(  0.0, 0.0,  20.0), Vector::new(0.0,  40.0,   0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Left:  back-bottom
    Triangle::new(Vector::new(10.0, 40.0, 100.0), Vector::new(  0.0, 0.0, -20.0), Vector::new(0.0, -40.0,   0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Left:  front-top
    Triangle::new(Vector::new(10.0,  0.0, 100.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0,  40.0,   0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Front: bottom-left
    Triangle::new(Vector::new(30.0, 40.0, 100.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0, -40.0,   0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Front: top-right
    Triangle::new(Vector::new(30.0,  0.0,  80.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0,  40.0,   0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Back:  bottom-right
    Triangle::new(Vector::new(10.0, 40.0,  80.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0, -40.0,   0.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Back:  top-left
    Triangle::new(Vector::new(10.0, 40.0, 100.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0,   0.0, -20.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Top:   front-left
    Triangle::new(Vector::new(30.0, 40.0,  80.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0,   0.0,  20.0), Color::new(0.0, 0.0, 0.0), Color::new(0.75, 0.75, 0.75)), // Top:   back-right
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

  let mut form_factors: HashMap<_, _> = (0..n).into_par_iter().map(|i| {
    let mut maps: Vec<HashMap<usize, Vec<f64>>> = vec![HashMap::new(); tris[i].patches.len()];

    for p in 0..tris[i].patches.len() {
      for j in 0..n {
        maps[p].insert(j, vec![0.0; tris[j].patches.len()]);
      }
    }

    (i, maps)
  }).collect();

  // Loop over all triangles in scene.
  for i in 0..n {
    print!("{} ", i);

    // Loop over all patches in rectangle i.
    for p_i in 0..tris[i].patches.len() {
      if p_i % tris[i].divisions as usize == 0 {
        print!("*");
        io::stdout().flush().unwrap();
      }

      // Loop over all triangles in scene for triangles i.
      for j in (i + 1)..n {
        // Loop over all patches in rectangle j.
        for p_j in 0..tris[j].patches.len() {
          // Monte Carlo integration of form factor double integral.

          let t_i: &PatchTriangle = &tris[i].sub_triangles[p_i];
          let t_j: &PatchTriangle = &tris[j].sub_triangles[p_j];

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
            let mut normal = Vector::new(0.0, 0.0, 0.0);
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

          form_factors.get_mut(&i).unwrap()[p_i].get_mut(&j).unwrap().insert(p_j, form_factor);
        }
      }
    }

    println!();
  }

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
          let area = tris[i].sub_triangles[p_i].area;
          let form_factor = form_factors[&i][p_i][&j][p_j] / area;
          form_factors.get_mut(&i).unwrap()[p_i].get_mut(&j).unwrap()[p_j] = if form_factor < 0.0 { 0.0 } else if form_factor > 1.0 { 1.0 } else { form_factor};
        }
      }
    }
  }

  form_factors
}

fn calculate_radiosity(tris: &mut [Triangle], form_factors: &HashMap<usize, Vec<HashMap<usize, Vec<f64>>>>) {
  for i in 0..tris.len() {
    for p_a in 0..tris[i].patches.len() {
      let color = (0..tris.len()).into_par_iter().map(|j| {
        if i == j { return Color::new(0.0, 0.0, 0.0); }
        (0..tris[j].patches.len()).into_iter().map(|p_b| {
          form_factors[&i][p_a][&j][p_b] * tris[j].patches[p_b]
        }).sum()
      }).sum();

      // Multiply sum with color of patch and add emission and
      // store overall patch radiosity of current iteration.
      tris[i].patches[p_a] = tris[i].color.entrywise_product(&color) + tris[i].emission;
    }
  }
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
      vertex_colors.get_mut(&tri.normal).unwrap().insert(tri.sub_triangles[p].a, Color::new(0.0, 0.0, 0.0));
      vertex_colors.get_mut(&tri.normal).unwrap().insert(tri.sub_triangles[p].b, Color::new(0.0, 0.0, 0.0));
      vertex_colors.get_mut(&tri.normal).unwrap().insert(tri.sub_triangles[p].c, Color::new(0.0, 0.0, 0.0));
      vertex_counts.get_mut(&tri.normal).unwrap().insert(tri.sub_triangles[p].a, 0);
      vertex_counts.get_mut(&tri.normal).unwrap().insert(tri.sub_triangles[p].b, 0);
      vertex_counts.get_mut(&tri.normal).unwrap().insert(tri.sub_triangles[p].c, 0);
    }
  }

  for tri in tris {
    for p in 0..tri.patches.len() {
      *vertex_counts.get_mut(&tri.normal).unwrap().get_mut(&tri.sub_triangles[p].a).unwrap() += 1;
      *vertex_counts.get_mut(&tri.normal).unwrap().get_mut(&tri.sub_triangles[p].b).unwrap() += 1;
      *vertex_counts.get_mut(&tri.normal).unwrap().get_mut(&tri.sub_triangles[p].c).unwrap() += 1;
    }
  }

  for tri in tris {
    for p in 0..tri.patches.len() {
      *vertex_colors.get_mut(&tri.normal).unwrap().get_mut(&tri.sub_triangles[p].a).unwrap() += tri.patches[p];
      *vertex_colors.get_mut(&tri.normal).unwrap().get_mut(&tri.sub_triangles[p].b).unwrap() += tri.patches[p];
      *vertex_colors.get_mut(&tri.normal).unwrap().get_mut(&tri.sub_triangles[p].c).unwrap() += tri.patches[p];
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
  let background_color = Color::new(0.0, 0.0, 0.0);

  let mut t = -1.0;
  let mut id = -1;
  let mut normal = Vector::new(0.0, 0.0, 0.0);

  // Find intersected triangle.
  if !intersect_scene(tris, ray, &mut t, &mut id, &mut normal) {
    return (background_color, background_color);
  }

  // Find intersected patch.
  let obj: &Triangle = &tris[id as usize];

  let find_patch_index = |triangle: &Triangle, ray: &Ray| -> Option<usize> {
    for p in 0..triangle.sub_triangles.len() {
      if triangle.sub_triangles[p].intersect(ray) != 0.0 {
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

  let bary = obj.sub_triangles[idx].barycentric_coordinates_at(&hitpoint);

  let a: Color = vertex_colors[&obj.normal][&obj.sub_triangles[idx].a];
  let b: Color = vertex_colors[&obj.normal][&obj.sub_triangles[idx].b];
  let c: Color = vertex_colors[&obj.normal][&obj.sub_triangles[idx].c];

  let interpolated: Color = a * bary.x.into_inner() + b * bary.z.into_inner() + c * bary.y.into_inner();

  return (obj.patches[idx], interpolated);
}


fn main() {
  let mut tris = tris();

  let divisions: u64 = 10;
  let samples = 10;

  println!("Calculating form factors ...");
  let form_factors = calculate_form_factors(&mut tris, divisions, samples);

  println!("Calculating radiosity ...");
  let iterations = 40;
  for i in 0..iterations {
    print!("{} ", i);
    calculate_radiosity(&mut tris, &form_factors);
  }
  println!();

  println!("Calculating vertex colors ...");
  let vertex_colors = calculate_vertex_colors(&tris);

  println!("Rendering images ...");

  let width = 640;
  let height = 480;

  let camera = Ray::new(&Vector::new(50.0, 52.0, 295.6), &Vector::new(0.0, -0.042612, -1.0).normalize());

  // Image edge vectors for pixel sampling.
  let cx: Vector = Vector::new(width as f64 * 0.5135 / height as f64, 0.0, 0.0);
  let cy: Vector = cx.cross_product(&camera.dir).normalize() * 0.5135;

  let mut image = Image::new(width, height);
  let mut image_interpolated = Image::new(width, height);

  let (tx, rx) = mpsc::channel();
  let sender = Mutex::new(tx);

  let start = Instant::now();

  // Loop over image rows.
  (0..height).into_par_iter().for_each(move |y| {
    let between = Range::new(0.0, 1.0);
    let mut rng = rand::thread_rng();

    // Loop over image columns.
    for x in 0..width {
      let mut accumulated_radiance = Color::new(0.0, 0.0, 0.0);
      let mut accumulated_radiance_interpolated = Color::new(0.0, 0.0, 0.0);

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
                              cy * (((y as f64) + ((sy as f64) + 0.5 + dy) / 2.0) / (height as f64) - 0.5) +
                              camera.dir;

            // Extend camera ray to start inside box.
            let start: Vector = camera.org + dir * 130.0;

            let (color, color_interpolated) = radiance(&tris, &Ray::new(&start, &dir.normalize()), &vertex_colors);

            // Determine constant radiance.
            accumulated_radiance += color / samples as f64;

            // Determine interpolated radiance.
            accumulated_radiance_interpolated += color_interpolated / samples as f64;
          }
        }
      }

      sender.lock().unwrap()
        .send((x, y, accumulated_radiance, accumulated_radiance_interpolated)).unwrap();
    }
  });

  for (x, y, radiance, radiance_interpolated) in rx {
    image.add_color(x, y, &radiance);
    image_interpolated.add_color(x, y, &radiance_interpolated);
  }

  let elapsed = start.elapsed();
  println!("Rendering took {} ms",
           (elapsed.as_secs() * 1_000) + (elapsed.subsec_nanos() / 1_000_000) as u64);


  println!("Saving images ...");

  if let Err(e) = image.save(&"image_patches.ppm".to_string()) {
    panic!("{:?}", e)
  }

  if let Err(e) = image_interpolated.save(&"image_smooth.ppm".to_string()) {
    panic!("{:?}", e)
  }
}
