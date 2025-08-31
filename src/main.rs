use physical_constants::{NEWTONIAN_CONSTANT_OF_GRAVITATION, SPEED_OF_LIGHT_IN_VACUUM};
use raylib::prelude::*;

pub const SCREEN_WIDTH: i32 = 800;
pub const SCREEN_HEIGHT: i32 = 600;

pub const G: f64 = NEWTONIAN_CONSTANT_OF_GRAVITATION;
pub const C: f64 = SPEED_OF_LIGHT_IN_VACUUM; //m/s
pub const SOLAR_MASS: f64 = 1.988416e30f64; //kg
pub const SCALE_FACTOR: f64 = 3.4e-11f64; //scale down
pub const SPEEDUP_FACTOR: f64 = 0.000001; // acceleration of frame so light movement is visible

struct Beem {
    position: Vector2,
    angle: f32, // 0 - 2PI
    trail: Vec<Vector2>,
    is_alive: bool,
    original_position: Vector2,
    original_angle: f32,

    r: f64,    // radial distance from black hole
    phi: f64,  // angular coordinate
    dr: f64,   // radial velocity
    dphi: f64, // angular velocity

    initialized: bool,
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

            r: 0.0,
            phi: 0.0,
            dr: 0.0,
            dphi: 0.0,
            initialized: false,
        }
    }

    fn initialize_geodesics(&mut self, black_hole: &BlackHole) {
        if self.initialized {
            return;
        }

        let dx = self.position.x - black_hole.position.x;
        let dy = self.position.y - black_hole.position.y;
        self.r = ((dx * dx + dy * dy) as f64).sqrt() / SCALE_FACTOR;
        self.phi = (dy as f64).atan2(dx as f64);

        let angle_to_bh = (dy as f64).atan2(dx as f64);
        let relative_angle = self.angle as f64 - angle_to_bh;

        self.dr = (C * relative_angle.cos()) / SCALE_FACTOR;
        self.dphi = (C * relative_angle.sin()) / (self.r * SCALE_FACTOR);

        self.initialized = true;
    }

    fn step(&mut self, black_hole: &BlackHole, r_s: f64, dt: f64) {
        if self.r < r_s * 1.01 {
            self.is_alive = false;
            return;
        }

        let d2r = self.r * self.dphi * self.dphi - (C * C * r_s) / (2.0 * self.r * self.r);
        let d2phi = -2.0 * self.dr * self.dphi / self.r;

        self.dr += d2r * dt;
        self.dphi += d2phi * dt;
        self.r += self.dr * dt;
        self.phi += self.dphi * dt;

        let x = (self.r * self.phi.cos()) * SCALE_FACTOR + black_hole.position.x as f64;
        let y = (self.r * self.phi.sin()) * SCALE_FACTOR + black_hole.position.y as f64;

        self.position = Vector2::new(x as f32, y as f32);
    }

    fn update(&mut self, d: &mut RaylibDrawHandle, black_hole: &BlackHole) {
        if !self.is_alive {
            if !self.trail.is_empty() {
                self.trail.remove(0);
            }
            return;
        }

        if !self.initialized {
            self.initialize_geodesics(black_hole);
        }

        self.trail.push(self.position);

        if self.trail.len() > Self::MAX_TRAIL_LENGTH {
            self.trail.remove(0);
        }

        let r_s_physical = black_hole.event_horizon_radius() / SCALE_FACTOR;
        if self.r < r_s_physical {
            self.is_alive = false;
            return;
        }

        let dt = d.get_frame_time() as f64 * SPEEDUP_FACTOR;
        self.step(black_hole, r_s_physical, dt);

        if self.position.x > SCREEN_WIDTH as f32
            || self.position.x < 0.0
            || self.position.y > SCREEN_HEIGHT as f32
            || self.position.y < 0.0
        {
            self.is_alive = false;
        }
    }

    pub fn should_remove(&self) -> bool {
        !self.is_alive && self.trail.is_empty()
    }

    pub fn respawn(&mut self) {
        self.position = self.original_position;
        self.angle = self.original_angle;
        self.trail.clear();
        self.is_alive = true;

        // Reset geodesic variables
        self.r = 0.0;
        self.phi = 0.0;
        self.dr = 0.0;
        self.dphi = 0.0;
        self.initialized = false;
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
}

impl Drop for Beem {
    fn drop(&mut self) {
        // destructor
    }
}

struct BlackHole {
    position: Vector2,
    mass: f64,   // in solar masses
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
        .title("Gargantua")
        .build();
    rl.set_target_fps(60);
    let sh: f32 = SCREEN_HEIGHT as f32;
    let sw: f32 = SCREEN_WIDTH as f32;

    // black hole definition, simulating Interstellar's black hole Gargantua
    let gargantua = BlackHole::new(sw / 2.0, sh / 2.0, 10e8);

    // beams creation logic with uniform spacing for visualization
    let mut beems: Vec<Beem> = Vec::new();
    let num_beams = 10;
    let beam_spacing = sh / (num_beams + 1) as f32;
    for i in 1..=num_beams {
        let y_pos = i as f32 * beam_spacing;
        beems.push(Beem::new(0.0, y_pos, 0.0));
    }

    // main loop that runs simulation until window closes
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);

        // draw background
        d.clear_background(Color::LIGHTGRAY);

        // update and draw beams, respawn when trail is completely gone
        for beem in &mut beems {
            beem.update(&mut d, &gargantua);
            beem.draw(&mut d);

            // respawn beam if it should be removed
            if beem.should_remove() {
                beem.respawn();
            }
        }

        // draw event horizon
        gargantua.draw(&mut d);

        // vars for description
        let mass_text = format!("mass: {:.1e} solar masses", gargantua.mass);
        let radius_km = gargantua.event_horizon_radius() / SCALE_FACTOR;
        let radius_text = format!("event horizon: {:.1} km", radius_km);
        let speedup_text = format!("speedup: x{:.1e}", SPEEDUP_FACTOR);

        // description with relevant info
        d.draw_text(&mass_text, 10, 10, 20, Color::BLACK);
        d.draw_text(&radius_text, 10, 35, 20, Color::BLACK);
        d.draw_text(&speedup_text, 10, 60, 20, Color::BLACK);
    }
}
