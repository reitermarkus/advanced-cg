extern crate rand;
extern crate crossbeam;
extern crate rayon;

#[macro_use]
extern crate lazy_static;

mod color;
mod vector;
mod triangle;
mod ray;
mod scene_object;
mod sphere;
mod image;
mod sampler;

use ray::Ray;
use triangle::Triangle;
use vector::Vector;
use color::Color;
use scene_object::SceneObject;
use scene_object::ReflType;
use sphere::Sphere;
use image::Image;

use sampler::{drand48, random_direction};
use rayon::prelude::*;

use std::f64::consts::PI;
use std::sync::Arc;
use std::thread;

lazy_static! {
  static ref TRIS: Vec<Triangle> = {
    vec![
      /* Cornell Box walls */
      Triangle::new(Vector::new(  0.0,  0.0,   0.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Back:   bottom-left
      Triangle::new(Vector::new(100.0, 80.0,   0.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Back:   top-right
      Triangle::new(Vector::new(  0.0,  0.0, 170.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,   0.0, -170.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Bottom: front-left
      Triangle::new(Vector::new(100.0,  0.0,   0.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,   0.0,  170.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Bottom: back-right
      Triangle::new(Vector::new(  0.0, 80.0,   0.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,   0.0,  170.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Top:    back-left
      Triangle::new(Vector::new(100.0, 80.0, 170.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,   0.0, -170.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Top:    front-right
      Triangle::new(Vector::new(  0.0,  0.0, 170.0), Vector::new(   0.0, 0.0, -170.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.75, 0.25, 0.25), ReflType::DIFF), // Left:   front-bottom
      Triangle::new(Vector::new(  0.0, 80.0,   0.0), Vector::new(   0.0, 0.0,  170.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.75, 0.25, 0.25), ReflType::DIFF), // Left:   back-top
      Triangle::new(Vector::new(100.0,  0.0,   0.0), Vector::new(   0.0, 0.0,  170.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.25, 0.25, 0.75), ReflType::DIFF), // Right:  back-bottom
      Triangle::new(Vector::new(100.0, 80.0, 170.0), Vector::new(   0.0, 0.0, -170.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.25, 0.25, 0.75), ReflType::DIFF), // Right:  front-top
      Triangle::new(Vector::new(100.0,  0.0, 170.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.25, 0.75,  0.25), ReflType::DIFF),  // Front:  bottom-right (not visible)
      Triangle::new(Vector::new(  0.0, 80.0, 170.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.25, 0.75,  0.25), ReflType::DIFF),  // Front:  top-left (not visible)

      /* Cuboid in room */
      Triangle::new(Vector::new(30.0,  0.0, 100.0), Vector::new(  0.0, 0.0, -20.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Right: front-bottom
      Triangle::new(Vector::new(30.0, 40.0,  80.0), Vector::new(  0.0, 0.0,  20.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Right: back-top
      Triangle::new(Vector::new(10.0,  0.0,  80.0), Vector::new(  0.0, 0.0,  20.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Left:  back-bottom
      Triangle::new(Vector::new(10.0, 40.0, 100.0), Vector::new(  0.0, 0.0, -20.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Left:  front-top
      Triangle::new(Vector::new(10.0,  0.0, 100.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Front: bottom-left
      Triangle::new(Vector::new(30.0, 40.0, 100.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Front: top-right
      Triangle::new(Vector::new(30.0,  0.0,  80.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Back:  bottom-right
      Triangle::new(Vector::new(10.0, 40.0,  80.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Back:  top-left
      Triangle::new(Vector::new(10.0, 40.0, 100.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0,   0.0, -20.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Top:   front-left
      Triangle::new(Vector::new(30.0, 40.0,  80.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0,   0.0,  20.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Top:   back-right
    ]
  };

  static ref SPHERES: Vec<Sphere> = {
    vec![
      Sphere::new(16.5, Vector::new(27.0, 16.5, 47.0), Color::zero(), Color::new(1.0, 1.0, 1.0), ReflType::REFR),  /* Mirror sphere */
      Sphere::new(16.5, Vector::new(73.0, 16.5, 78.0), Color::zero(), Color::new(1.0, 1.0, 1.0), ReflType::SPEC),  /* Mirror sphere */
      Sphere::new(1.5, Vector::new(50.0, 81.6 - 16.5, 81.6), Color::new(4.0, 4.0, 4.0) * 100.0, Color::zero(), ReflType::DIFF) /* Light */
    ]
  };
}

fn intersect_scene(objects: &[Box<&SceneObject>], ray: &Ray, t: &mut f64, id: &mut usize) -> bool {
  *t = 1e20;

  for i in 0..objects.len() {
    let d = objects[i].intersect(ray);
    if d > 0.0 && d < *t {
      *t = d;
      *id = i;
    }
  }

  *t < 1e20
}

fn perfect_reflection(dir: &Vector, normal: &Vector) -> Vector {
  *dir - *normal * 2.0 * (*normal).dot_product(dir)
}

fn radiance(objects: &[Box<&SceneObject>], ray: &Ray, mut depth: i32, E: i32) -> Color {
  depth += 1;
  let mut t = 0.0;
  let mut id = 0;

  if !intersect_scene(objects, ray, &mut t, &mut id) {
     /* No intersection with scene */
    return Color::new(0.0, 0.0, 0.0);
  }

  /* Intersection point */
  let hitpoint = ray.org + ray.dir + t;
  let obj = &objects[id];

  /* Normal at intersection */
  let mut normal: Vector;

  if obj.is_sphere() {
    let sphere: &Sphere = match obj.as_any().downcast_ref::<Sphere>() {
      Some(s) => s,
      None => panic!("obj is not a sphere!"),
    };

    normal = (hitpoint - sphere.position).normalize();
  } else {
    let triangle: &Triangle = match obj.as_any().downcast_ref::<Triangle>() {
      Some(tr) => tr,
      None => panic!("obj is not a triangle!"),
    };

    normal = triangle.normal;
  }

  let mut nl = normal;

  /* Obtain flipped normal, if object hit from inside */
  if normal.dot_product(&ray.dir) >= 0.0 {
    nl = -nl;
  }

  let mut col = obj.color();

  /* Maximum RGB reflectivity for Russian Roulette */
  let p = col.clamp(0.0, 0.999).max();

  /* After 5 bounces or if max reflectivity is zero */
  if depth > 5 || p == 0.0 {
    if drand48(0.0, 1.0) >= p {
      /* No further bounces, only return potential emission */
      return obj.emission() * E as f64;
    }

    /* Scale estimator to remain unbiased */
    col = col * (1.0 / p);
  }

  let mut ray_dir = ray.dir;
  let mut nt = 1.5;

  match obj.refl() {
    ReflType::DIFF => {
      /* Compute random reflection vector on hemisphere */
      let r1 = 2.0 * PI * drand48(0.0, 1.0);
      let r2 = drand48(0.0, 1.0);
      let r2s = r2.sqrt();

      /* Set up local orthogonal coordinate system u,v,w on surface */
      let w = nl;
      let u = if w.x.abs() > 0.1 { Vector::new(0.0, 1.0, 0.0) } else { Vector::new(1.0, 0.0, 0.0) };
      let u = u.cross_product(&w).normalize();

      let v = w.cross_product(&u);

      /* Random reflection vector d */
      let d = (u * r1.cos() * r2s +
               v * r1.sin() * r2s +
               w * (1.0 - r2).sqrt()).normalize();

      let mut e: Vector = Vector::zero();
      for i in 0..objects.len() {
        let light_source = &objects[i];

        if light_source.emission().x <= 0.0 && light_source.emission().y <= 0.0 && light_source.emission().z <= 0.0 {
          continue;
        }

        if !light_source.is_sphere() {
          println!("Warning: Only spherical light sources are implemented.");
          continue;
        }

        let sphere: &Sphere = match light_source.as_any().downcast_ref::<Sphere>() {
          Some(s) => s,
          None => panic!("obj is not a sphere!"),
        };

        /* Randomly sample spherical light source from surface intersection */

        // Create random sample direction towards spherical light source.
        let cos_a_max = (1.0 - sphere.radius.powf(2.0) /
                                (hitpoint - sphere.position).dot_product(&(hitpoint - sphere.position))).sqrt();

        let l = random_direction(sphere.position - hitpoint, (cos_a_max).acos());

        if intersect_scene(objects, &Ray::new(hitpoint, l), &mut t, &mut id) && id == i {
          let omega = 2.0 * PI * (1.0 - cos_a_max);
          e += col.entrywise_product(sphere.emission() * l.dot_product(&nl) * omega) / PI;
        }
      }

          /* Return potential light emission, direct lighting, and indirect lighting (via
        recursive call for Monte-Carlo integration */
      return obj.emission() * E as f64 + e + col.entrywise_product(radiance(objects, &Ray::new(hitpoint, d), depth, 0));
    },
    ReflType::GLOS => {
      let angle = 0.10;
      let l = random_direction(nl, angle);

      return obj.emission() + col.entrywise_product(radiance(objects, &Ray::new(hitpoint, l), depth, 1));
    },
    ReflType::SPEC => {
      // Return light emission mirror reflection (via recursive call using perfect reflection vector).
      let reflection_ray = Ray::new(hitpoint, perfect_reflection(&ray.dir, &normal));
      return obj.emission() + col.entrywise_product(radiance(objects, &reflection_ray, depth, 1));
    },
    ReflType::TRAN => {
      let angle = 0.2;
      ray_dir = random_direction(ray.dir, angle);
      nt = 1.15;
    },
    _ => {}
  }

  let nc = 1.0;

  let into = normal.dot_product(&nl) > 0.0;       /* Bool for checking if ray from outside going in */
  let nnt = if into { nc / nt } else { nt / nc }; /* Set ratio depending on hit from inside or outside */

  let ddn = ray_dir.dot_product(&nl);
  let cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

  let reflection_ray = Ray::new(hitpoint, perfect_reflection(&ray_dir, &normal)); // Perfect reflection.

    // Check for total internal reflection, if so only reflect.
  if cos2t < 0.0 {
    return obj.emission() +  col.entrywise_product(radiance(objects, &reflection_ray, depth, 1));
  }


  /* Otherwise reflection and/or refraction occurs */
  let mut transmission_direction;

  // Determine transmitted ray direction for refraction.
  if into {
    transmission_direction = ray_dir * nnt - normal * (ddn * nnt + cos2t.sqrt());
  } else {
    transmission_direction = ray_dir * nnt + normal * (ddn * nnt + cos2t.sqrt());
  }

  let transmission_direction = transmission_direction.normalize();

  /* Cosine of correct angle depending on outside/inside */
  let c = if into { 1.0 + ddn } else { 1.0 - transmission_direction.dot_product(&normal) };

  /* Compute Schlick's approximation of Fresnel equation */
  let R0 = ((nc - nt) / (nc + nt)).powf(2.0);
  let Re = R0 + (1.0 - R0) * c.powf(5.0);   // Reflectance
  let Tr = 1.0 - Re;                     // Transmittance

  let transmission_ray = Ray::new(hitpoint, transmission_direction);

  // Initially, use both reflection and trasmission.
  if depth < 3 {
    return obj.emission() +  col.entrywise_product(radiance(objects, &reflection_ray, depth, 1) * Re
                                                   + radiance(objects, &transmission_ray, depth, 1) * Tr);
  }
  // Probability for selecting reflectance or transmittance */
  let P = 0.25 + 0.5 * Re;
  let RP = Re / P;         /* Scaling factors for unbiased estimator */
  let TP = Tr / (1.0 - P);

  /* Russian Roulette */
  if drand48(0.0, 1.0) < P {
    return obj.emission() +  col.entrywise_product(radiance(objects, &reflection_ray, depth, 1) * RP);
  }

  return obj.emission() +  col.entrywise_product(radiance(objects, &transmission_ray, depth, 1) * TP);
}

fn main() {
  let width = 1024;
  let height = 768;
  let samples = 2;

  let aperture = 2.6;
  let focal_length = 120.0;

  let mut scene_objects: Vec<Box<&SceneObject>> = Vec::new();

  for tris in TRIS.iter() {
    scene_objects.push(Box::new(tris));
  }

  for sphere in SPHERES.iter() {
    scene_objects.push(Box::new(sphere));
  }

  let camera = Ray::new(Vector::new(50.0, 52.0, 295.6), Vector::new(0.0, -0.042612, -1.0).normalize());
  let focal_point = camera.org + camera.dir * focal_length;

  let cx = Vector::new(width as f64 * 0.5135 / height as f64, 0.0, 0.0);
  let cy = (cx.cross_product(&camera.dir)).normalize() * 0.5135;

  let image = Arc::new(Image::new(width, height));
  let image_clone = image.clone();

  let image_thread = thread::spawn(move || {
    thread::park();

    if let Err(e) = image_clone.save("image.ppm") {
      panic!("{:?}", e)
    }
  });

  let start_saving = if height < 20 { height / 2 } else { height - 20 };

  // Loop over image rows.
  (0..height).into_par_iter().for_each(|y| {
     println!("\rRendering ({}spp) {}%     ", samples * 4, (100 * y / (height - 1)));

    // Loop over image columns.
    for x in 0..width {
      image.set_color(x, y, Color::zero());
      let mut accumulated_radiance = Color::zero();

      // 2 x 2 subsampling per pixel.
      for sy in 0..2 {
        for sx in 0..2 {

          // Computes radiance at subpixel using multiple samples.
          for _ in 0..samples {
            // Generate random sample on circular lens.
            let random_radius = drand48(0.0, 1.0);
            let random_angle = drand48(0.0, 1.0);
            let lens_sample_point = aperture * Vector::new(random_radius.sqrt() * (2.0 * PI * random_angle).cos(), random_radius.sqrt() * (2.0 * PI * random_angle).sin(), 0.0);

            let mut dir = (focal_point - (camera.org + lens_sample_point)).normalize();

            dir = (camera.dir + dir).normalize();

            let mut nu_filter_samples = || -> f64 {
              let r = 2.0 * drand48(0.0, 1.0) as f64;
              if r < 1.0 { r.sqrt() - 1.0 } else { 1.0 - (2.0 - r).sqrt() }
            };

            let dx = nu_filter_samples();
            let dy = nu_filter_samples();

            // Ray direction into scene from camera through sample.
            dir = cx * ((x as f64 + (sx as f64 + 0.5 + dx) / 2.0) / width as f64 - 0.5) +
                              cy * ((y as f64 + (sy as f64 + 0.5 + dy) / 2.0) / height as f64 - 0.5) +
                              dir;

            // Extend camera ray to start inside box.
            let start: Vector = camera.org + dir * 130.0;

            dir = dir.normalize();

            let ray = Ray::new(start + lens_sample_point, dir);

            /* Accumulate radiance */
            accumulated_radiance = accumulated_radiance + radiance(&scene_objects, &ray, 0, 1) / samples as f64;
          }

          accumulated_radiance = accumulated_radiance.clamp(0.0, 1.0) * 0.25;
        }
      }

      image.set_color(x, y, accumulated_radiance / samples as f64);

      if y == start_saving {
        image_thread.thread().unpark();
      }
    }
  });

  image_thread.join().unwrap();
}
