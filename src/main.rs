use physical_constants::{NEWTONIAN_CONSTANT_OF_GRAVITATION, SPEED_OF_LIGHT_IN_VACUUM};
use raylib::prelude::*;

pub const SCREEN_WIDTH: i32 = 800;
pub const SCREEN_HEIGHT: i32 = 600;

pub const G: f64 = NEWTONIAN_CONSTANT_OF_GRAVITATION;
pub const C: f64 = SPEED_OF_LIGHT_IN_VACUUM; //m/s
pub const SOLAR_MASS: f64 = 1.988416e30f64; //kg
pub const SCALE_FACTOR: f64 = 3.4e-11f64; //scale down
pub const SPEEDUP_FACTOR: f64 = 39235.0 / 2.0; // acceleration of frame so light movement is visible -> if 1, real time

struct Beem {
    position: Vector2,
    angle: f32, // 0 - 2PI
    trail: Vec<Vector2>,
    is_alive: bool,
    original_position: Vector2,
    original_angle: f32,
}

impl Beem {
    const RADIUS: f32 = 2.0;
    const MAX_TRAIL_LENGTH: usize = 100;

    pub fn new(x: f32, y: f32, angle: f32) -> Self {
        Self {
            position: Vector2::new(x, y),
            angle,
            trail: Vec::new(),
            is_alive: true,
            original_position: Vector2::new(x, y),
            original_angle: angle,
        }
    }

    fn update(&mut self, d: &mut RaylibDrawHandle, black_hole: &BlackHole) {
        if !self.is_alive {
            // If beam is dead, gradually remove trail points to make it fade away
            if !self.trail.is_empty() {
                self.trail.remove(0);
            }
            return;
        }

        // add current position to trail before moving
        self.trail.push(self.position);

        // limit trail length
        if self.trail.len() > Self::MAX_TRAIL_LENGTH {
            self.trail.remove(0);
        }

        // check if within event horizon before moving
        if self.within_event_horizon(black_hole) {
            self.is_alive = false;
            return;
        }

        let distance: f32 = (d.get_frame_time() as f64 * C * SCALE_FACTOR * SPEEDUP_FACTOR) as f32;
        self.position.x += distance * self.angle.cos();
        self.position.y += distance * self.angle.sin();

        // Instead of wrapping, kill the beam when it goes off screen
        if self.position.x > SCREEN_WIDTH as f32
            || self.position.x < 0.0
            || self.position.y > SCREEN_HEIGHT as f32
            || self.position.y < 0.0
        {
            self.is_alive = false;
            return;
        }
    }

    //collision detection
    fn within_event_horizon(&self, black_hole: &BlackHole) -> bool {
        let distance = ((self.position.x - black_hole.position.x).powi(2)
            + (self.position.y - black_hole.position.y).powi(2))
        .sqrt();
        distance <= black_hole.event_horizon_radius() as f32
    }

    pub fn draw(&self, d: &mut RaylibDrawHandle) {
        // draw trail
        for i in 1..self.trail.len() {
            let alpha = (i as f32 / self.trail.len() as f32 * 255.0) as u8;
            let trail_color = Color::new(255, 255, 255, alpha);
            d.draw_line_ex(
                self.trail[i - 1],
                self.trail[i],
                Self::RADIUS * 2.0,
                trail_color,
            );
        }

        // draw current position only if beam is still alive
        if self.is_alive {
            d.draw_circle_v(self.position, Self::RADIUS, Color::WHITE);
        }
    }

    // check if beam should be removed (no trail left and dead)
    pub fn should_remove(&self) -> bool {
        !self.is_alive && self.trail.is_empty()
    }

    // respawn the beam at original position
    pub fn respawn(&mut self) {
        self.position = self.original_position;
        self.angle = self.original_angle;
        self.trail.clear();
        self.is_alive = true;
    }
}

impl Drop for Beem {
    fn drop(&mut self) {
        // destructor
    }
}

struct BlackHole {
    position: Vector2,
    mass: f64,   //in solar masses
    radius: f64, // in pixels
}

impl BlackHole {
    pub fn new(x: f32, y: f32, mass: f64) -> Self {
        Self {
            position: Vector2::new(x, y),
            mass,
            radius: Self::event_horizon_radius_init(mass),
        }
    }
    // public r_s
    pub fn event_horizon_radius(&self) -> f64 {
        let mass_kg: f64 = self.mass * SOLAR_MASS;
        let r_s: f64 = SCALE_FACTOR * (2.0 * G * mass_kg) / (C * C);
        r_s
    }

    // r_s used for init
    fn event_horizon_radius_init(mass: f64) -> f64 {
        let mass_kg: f64 = mass * SOLAR_MASS;
        let r_s: f64 = SCALE_FACTOR * (2.0 * G * mass_kg) / (C * C);
        r_s
    }

    pub fn draw(&self, d: &mut RaylibDrawHandle) {
        let r_s = self.event_horizon_radius();
        d.draw_circle_v(self.position, r_s as f32, Color::BLACK);
    }
}

fn main() {
    let (mut rl, thread) = raylib::init()
        .size(SCREEN_WIDTH, SCREEN_HEIGHT)
        .title("Gargantula")
        .build();
    rl.set_target_fps(60);
    let sh: f32 = SCREEN_HEIGHT as f32;
    let sw: f32 = SCREEN_WIDTH as f32;

    // black hole definition, simmulating Interstrellar's black hole Gargantula
    let gargantula = BlackHole::new(sw / 2.0, sh / 2.0, 10e8);

    // beems creation logic with uniform spacing for visualization
    let mut beems: Vec<Beem> = Vec::new();
    let num_beams = 10;
    let beam_spacing = sh / (num_beams + 1) as f32;
    for i in 1..=num_beams {
        let y_pos = i as f32 * beam_spacing;
        beems.push(Beem::new(0.0, y_pos, 0.0));
    }

    // main loop that runs simmulation until window closes
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);

        // draw background
        d.clear_background(Color::LIGHTGRAY);

        //update and draw beems, respawn when trail is completely gone
        for beem in &mut beems {
            beem.update(&mut d, &gargantula);
            beem.draw(&mut d);

            // respawn beam if it should be removed
            if beem.should_remove() {
                beem.respawn();
            }
        }

        //draw event horizon
        gargantula.draw(&mut d);

        // vars for description
        let mass_text = format!("mass: {:.1e} solar masses", gargantula.mass);
        let radius_km = gargantula.event_horizon_radius() / SCALE_FACTOR;
        let radius_text = format!("event horizon: {:.1} km", radius_km);
        let speedup_text = format!("speedup: x{:.1e}", SPEEDUP_FACTOR);

        // description with relevant info
        d.draw_text(&mass_text, 10, 10, 20, Color::BLACK);
        d.draw_text(&radius_text, 10, 35, 20, Color::BLACK);
        d.draw_text(&speedup_text, 10, 60, 20, Color::BLACK);
    }
}
