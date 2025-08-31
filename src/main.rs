use physical_constants::{NEWTONIAN_CONSTANT_OF_GRAVITATION, SPEED_OF_LIGHT_IN_VACUUM};
use raylib::prelude::*;

pub const SCREEN_WIDTH: i32 = 800;
pub const SCREEN_HEIGHT: i32 = 600;
pub const NUM_BEEMS: i32 = 50;

pub const G: f64 = NEWTONIAN_CONSTANT_OF_GRAVITATION;
pub const C: f64 = SPEED_OF_LIGHT_IN_VACUUM; //m/s
pub const SOLAR_MASS: f64 = 1.988416e30f64; //kg
pub const SCALE_FACTOR: f64 = 3.4e-11f64; //scale down
pub const SPEEDUP_FACTOR: f64 = 0.0000015; // acceleration of frame so light movement is visible

#[derive(Copy, Clone)]
struct BeemState {
    r: f64,
    phi: f64,
    dr: f64,
    dphi: f64,
}

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

    fn get_derivatives(&self, state: BeemState, black_hole: &BlackHole) -> BeemState {
        let r_s = black_hole.event_horizon_radius() / SCALE_FACTOR;

        let d2r = state.r * state.dphi * state.dphi * (1.0 - r_s / state.r)
            - (C * C * r_s) / (2.0 * state.r * state.r);
        let d2phi = -2.0 * state.dr * state.dphi / state.r;

        BeemState {
            r: state.dr,
            phi: state.dphi,
            dr: d2r,
            dphi: d2phi,
        }
    }

    fn rk4_step(&self, state: BeemState, dt: f64, black_hole: &BlackHole) -> BeemState {
        let k1 = self.get_derivatives(state, black_hole);
        let k2_state = BeemState {
            r: state.r + k1.r * dt / 2.0,
            phi: state.phi + k1.phi * dt / 2.0,
            dr: state.dr + k1.dr * dt / 2.0,
            dphi: state.dphi + k1.dphi * dt / 2.0,
        };
        let k2 = self.get_derivatives(k2_state, black_hole);

        let k3_state = BeemState {
            r: state.r + k2.r * dt / 2.0,
            phi: state.phi + k2.phi * dt / 2.0,
            dr: state.dr + k2.dr * dt / 2.0,
            dphi: state.dphi + k2.dphi * dt / 2.0,
        };
        let k3 = self.get_derivatives(k3_state, black_hole);

        let k4_state = BeemState {
            r: state.r + k3.r * dt,
            phi: state.phi + k3.phi * dt,
            dr: state.dr + k3.dr * dt,
            dphi: state.dphi + k3.dphi * dt,
        };
        let k4 = self.get_derivatives(k4_state, black_hole);

        BeemState {
            r: state.r + (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r) * dt / 6.0,
            phi: state.phi + (k1.phi + 2.0 * k2.phi + 2.0 * k3.phi + k4.phi) * dt / 6.0,
            dr: state.dr + (k1.dr + 2.0 * k2.dr + 2.0 * k3.dr + k4.dr) * dt / 6.0,
            dphi: state.dphi + (k1.dphi + 2.0 * k2.dphi + 2.0 * k3.dphi + k4.dphi) * dt / 6.0,
        }
    }

    fn step(&mut self, black_hole: &BlackHole, r_s: f64, dt: f64) {
        if self.r < r_s * 1.01 {
            self.is_alive = false;
            return;
        }

        let current_state = BeemState {
            r: self.r,
            phi: self.phi,
            dr: self.dr,
            dphi: self.dphi,
        };

        let new_state = self.rk4_step(current_state, dt, black_hole);

        self.r = new_state.r;
        self.phi = new_state.phi;
        self.dr = new_state.dr;
        self.dphi = new_state.dphi;

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

        self.r = 0.0;
        self.phi = 0.0;
        self.dr = 0.0;
        self.dphi = 0.0;
        self.initialized = false;
    }

    pub fn draw(&self, d: &mut RaylibDrawHandle) {
        for i in 1..self.trail.len() {
            let alpha = (i as f32 / self.trail.len() as f32 * 255.0) as u8;
            let trail_color = Color::new(169, 169, 169, alpha);
            d.draw_line_ex(
                self.trail[i - 1],
                self.trail[i],
                Self::RADIUS * 2.0,
                trail_color,
            );
        }

        if self.is_alive {
            d.draw_circle_v(self.position, Self::RADIUS, Color::DARKGRAY);
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
    mass: f64,
    radius: f64,
}

impl BlackHole {
    pub fn new(x: f32, y: f32, mass: f64) -> Self {
        Self {
            position: Vector2::new(x, y),
            mass,
            radius: Self::event_horizon_radius_init(mass),
        }
    }

    pub fn event_horizon_radius(&self) -> f64 {
        let mass_kg: f64 = self.mass * SOLAR_MASS;
        let r_s: f64 = SCALE_FACTOR * (2.0 * G * mass_kg) / (C * C);
        r_s
    }

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

    let gargantua = BlackHole::new(sw / 2.0, sh / 2.0, 10e8);

    let mut beems: Vec<Beem> = Vec::new();
    let num_beams = NUM_BEEMS;
    let beam_spacing = sh / (num_beams + 1) as f32;
    for i in 1..=num_beams {
        let y_pos = i as f32 * beam_spacing;
        beems.push(Beem::new(0.0, y_pos, 0.0));
    }

    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);

        d.clear_background(Color::new(20, 20, 20, 255));

        for beem in &mut beems {
            beem.update(&mut d, &gargantua);
            beem.draw(&mut d);

            if beem.should_remove() {
                beem.respawn();
            }
        }

        gargantua.draw(&mut d);

        let mass_text = format!("mass: {:.1e} solar masses", gargantua.mass);
        let light_beems_text = format!("number of light beems: {:}", NUM_BEEMS);
        let radius_km = gargantua.event_horizon_radius() / SCALE_FACTOR;
        let radius_text = format!("event horizon: {:.1} km", radius_km);
        let speedup_text = format!("speedup: x{:.1e}", SPEEDUP_FACTOR);

        d.draw_text(&mass_text, 10, 10, 20, Color::WHITE);
        d.draw_text(&radius_text, 10, 35, 20, Color::WHITE);
        d.draw_text(&light_beems_text, 10, 60, 20, Color::WHITE);
        d.draw_text(&speedup_text, 10, 85, 20, Color::WHITE);
    }
}
