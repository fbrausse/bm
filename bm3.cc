
#include <cmath>
#include <cinttypes>
#include <utility>	/* std::pair */
#include <queue>	/* std::priority_queue */
#include <vector>	/* std::vector */
#include <memory>	/* std::shared_ptr */
#include <iostream>
#include <random>
#include <algorithm>	/* std::none_of */
#include <set>
#include <cstring>
#include <unistd.h>	/* getopt(3) */

template <typename T>
struct vec2 {
	T x, y;

	friend vec2 operator+(const vec2 &a, const vec2 &b) noexcept { return { a.x+b.x, a.y+b.y }; }
	friend vec2 operator-(const vec2 &a, const vec2 &b) noexcept { return a+-b; }
	friend vec2 operator-(const vec2 &a)                noexcept { return { -a.x, -a.y }; }
	friend vec2 operator*(const vec2 &a, const T &b)    noexcept { return { a.x * b, a.y * b }; }
	friend vec2 operator*(const T &a, const vec2 &b)    noexcept { return b*a; }
	friend vec2 operator/(const vec2 &a, const T &b)    noexcept { return { a.x / b, a.y / b }; }

	friend vec2 & operator+=(vec2 &a, const vec2 &b) noexcept { return a = a+b; }
	friend vec2 & operator-=(vec2 &a, const vec2 &b) noexcept { return a = a-b; }
	friend vec2 & operator*=(vec2 &a, const T &b)    noexcept { return a = a*b; }
	friend vec2 & operator/=(vec2 &a, const T &b)    noexcept { return a = a/b; }

	friend T dot(const vec2 &a, const vec2 &b) noexcept { return a.x*b.x + a.y*b.y; }
	friend T norm_sq(const vec2 &a)            noexcept { return dot(a,a); }
	friend T norm(const vec2 &a)               noexcept { using std::sqrt; return sqrt(norm_sq(a)); }
	friend T arg(const vec2 &a)                noexcept { using std::atan2; return atan2(a.y, a.x); }

	friend std::ostream & operator<<(std::ostream &os, const vec2 &a) noexcept
	{ return os << "[" << a.x << ";" << a.y << "]"; }
};

typedef vec2<double> vec2d;

struct Particle {
	vec2d p, v;
	double m, r;

	Particle() noexcept {}
	Particle(const vec2d &p, const vec2d &v, double m, double r)
	noexcept : p(p), v(v), m(m), r(r) {}

	void advance(double t) noexcept { p += v * t; }

	friend std::ostream & operator<<(std::ostream &os, const Particle &p) noexcept
	{
		return os << "Ptcl:p="<<p.p<<",v="<<p.v<<",m:"<<p.m<<",r:"<<p.r;
	}
};

static void hsv2rgb(double h, double s, double v, double *r, double *g, double *b)
{
	double f = fmod(h * 360.0, 60.0) / 60.0;
	double p = v * (1.0 - s);
	double q = v * (1.0 - f * s);
	double t = v * (1.0 - (1.0 - f) * s);
	switch ((int)(h * 6.0L)) {  // the long double (h * 6.0L) is important, otherwise h = 0.5 will result in 2 instead of 3
	case 0: *r = v; *g = t; *b = p; break;
	case 1: *r = q; *g = v; *b = p; break;
	case 2: *r = p; *g = v; *b = t; break;
	case 3: *r = p; *g = q; *b = v; break;
	case 4: *r = t; *g = p; *b = v; break;
	case 5: *r = v; *g = p; *b = q; break;
	}
}

struct ColoredParticle : public Particle {
	double red, green, blue;

	ColoredParticle() noexcept : Particle() { set_hsv(0.0, 0.0, 0.0); }

	ColoredParticle(const vec2d &p, const vec2d &v, double m, double r, double hue)
	noexcept : Particle(p, v, m, r) { set_hsv(hue); }

	void set_hsv(double h, double s = 1.0, double v = 1.0) noexcept
	{
		hsv2rgb(h, s, v, &red, &green, &blue);
	}
};

static bool collision_swap_colors = false;

static void collide(Particle &a, Particle &b) noexcept
{
	double m = a.m + b.m;
	vec2d p = a.p - b.p;
	vec2d v = a.v - b.v;
	double pv = dot(v,p);
	double pvn = pv / norm_sq(p);
	a.v -= 2*b.m/m * pvn * p;
	b.v += 2*a.m/m * pvn * p;
	if (collision_swap_colors) {
		ColoredParticle &ca = static_cast<ColoredParticle &>(a);
		ColoredParticle &cb = static_cast<ColoredParticle &>(b);
		std::swap(ca.red  , cb.red);
		std::swap(ca.green, cb.green);
		std::swap(ca.blue , cb.blue);
	}
}

static std::pair<bool,double> time_to_collision(const Particle &a, const Particle &b) noexcept
{
	using std::sqrt;

	vec2d  p   = b.p - a.p;
	double r   = a.r + b.r;
	double pr2 = norm_sq(p) - r*r;
	if (pr2 < 0)
		return { false, INFINITY };
	vec2d  v   = b.v - a.v;
	double nv2 = norm_sq(v);
	double pv  = dot(p,v);
	double D   = pv*pv - nv2*pr2;
	if (D < 0)
		return { false, INFINITY };
	double d   = D > 0 ? sqrt(D) : 0;
	double t1  = -pv - d;
	double t2  = -pv + d;
	double t   = true || t1 > 0 ? t1 : t2;
//	std::cerr << "collide " << a << " vs " << b << ": " << t1 << " " << t2 << std::endl;
	return { t > 0, t/nv2 };
}

struct Ev {
	double t;

	explicit Ev(double t) noexcept : t(t) {}
	virtual ~Ev() {}

	virtual void handle() const noexcept = 0;
	virtual bool handles(const Particle &) const noexcept = 0;
	virtual void invalidate(std::set<Particle *> &recompute) {}

	friend bool operator< (const Ev &a, const Ev &b) { return a.t < b.t; }
	friend bool operator<=(const Ev &a, const Ev &b) { return !(b < a); }
	friend bool operator> (const Ev &a, const Ev &b) { return   b < a ; }
	friend bool operator>=(const Ev &a, const Ev &b) { return !(a < b); }

	friend std::ostream & operator<<(std::ostream &os, const Ev &ev) noexcept
	{
		ev.print(os);
		return os;
	}

protected:
	virtual void print(std::ostream &os) const noexcept
	{
		os << typeid(*this).name() << " @ t = " << t;
	}
};

template <typename T>
struct shared_less {
	bool operator()(
		const std::shared_ptr<T> &a,
		const std::shared_ptr<T> &b
	) const noexcept {
		return *a < *b;
	}
};

using PQ = std::set<std::shared_ptr<Ev>,shared_less<Ev>>;

static bool overlaps(const Particle &a, const Particle &b)
{
	double r = a.r + b.r;
	return norm_sq(a.p - b.p) <= r*r;
}

struct Field {
	const double w, h;
	std::vector<ColoredParticle> particles;

	Field(double w, double h) : w(w), h(h) {}

	void advance(double t) noexcept { for (Particle &a : particles) a.advance(t); }

	bool fits(const Particle &p) const
	{
		return 0 <= p.p.x - p.r && p.p.x + p.r < w &&
		       0 <= p.p.y - p.r && p.p.y + p.r < h &&
		       std::none_of(particles.begin(), particles.end(),
		                    [&p](const Particle &q) { return overlaps(p,q); });
	}

	void next_collisions(PQ &events);
	void next_collisions(Particle &a, PQ &events);
	std::shared_ptr<Ev> next_collision(Particle &a, PQ &events);

	friend std::ostream & operator<<(std::ostream &os, const Field &f)
	{
		os << "Field:w=" << f.w << ";h=" << f.h << ";[\n";
		for (const Particle &a : f.particles)
			os << "\t" << a << "\n";
		return os << "]";
	}
};

struct C1Ev : public Ev {
	Particle &a;
	const bool vert;
	Field &f;
	PQ &events;

	C1Ev(double t, Particle &a, bool vert, Field &f, PQ &events) noexcept
	: Ev(t), a(a), vert(vert), f(f), events(events) {}

	bool handles(const Particle &p) const noexcept { return &a == &p; }

	void handle() const noexcept
	{
		if (vert) a.v.x = -a.v.x; else a.v.y = -a.v.y;

		std::set<Particle *> rem = { &a };
		for (auto it = events.begin(); it != events.end();)
			if ((*it)->handles(a)) {
				(*it)->invalidate(rem);
				it = events.erase(it);
			} else
				++it;
		for (Particle *q : rem) {
			std::shared_ptr<Ev> e = f.next_collision(*q, events);
			if (e)
				events.insert(e);
		}
	}

	void print(std::ostream &os) const noexcept
	{
		os << "C1Ev(" << t << ", " << a << ", " << vert << ")";
	}
};

struct C2Ev : public Ev {
	Particle &a, &b;
	Field &f;
	PQ &events;

	C2Ev(double t, Particle &a, Particle &b, Field &f, PQ &events) noexcept
	: Ev(t), a(a), b(b), f(f), events(events) {}

	bool handles(const Particle &p) const noexcept { return &a == &p || &b == &p; }

	void invalidate(std::set<Particle *> &recompute) { recompute.insert(&a); recompute.insert(&b); }

	void handle() const noexcept
	{
		collide(a, b);

		std::set<Particle *> rem = { &a, &b };
		for (auto it = events.begin(); it != events.end();)
			if ((*it)->handles(a) || (*it)->handles(b)) {
				(*it)->invalidate(rem);
				it = events.erase(it);
			} else
				++it;
		for (Particle *q : rem) {
			std::shared_ptr<Ev> e = f.next_collision(*q, events);
			if (e)
				events.insert(std::move(e));
		}
	}

	void print(std::ostream &os) const noexcept
	{
		os << "C2Ev(" << t << ", " << a << ", " << b << ")";
	}
};

void Field::next_collisions(PQ &events)
{
	for (Particle &a : particles) {
		const vec2d &v = a.v;
		double tx = v.x < 0 ? std::max(0.,  a.p.x-a.r)/-v.x
			  : v.x > 0 ? std::min(w ,w-a.p.x-a.r)/ v.x
			  : INFINITY;
		if (tx < INFINITY)
			events.insert(std::make_shared<C1Ev>(tx, a, true, *this, events));
		double ty = v.y < 0 ? std::max(0.,  a.p.y-a.r)/-v.y
			  : v.y > 0 ? std::min(h ,h-a.p.y-a.r)/ v.y
			  : INFINITY;
		if (ty < INFINITY)
			events.insert(std::make_shared<C1Ev>(ty, a, false, *this, events));
	}

	for (unsigned i=0; i<particles.size(); i++) {
		Particle &a = particles[i];
		for (unsigned j=i+1; j<particles.size(); j++) {
			Particle &b = particles[j];
			if (overlaps(a, b))
				continue;
			auto r = time_to_collision(a, b);
			if (!r.first)
				continue;
			events.insert(std::make_shared<C2Ev>(r.second, a, b, *this, events));
		}
	}
}

void Field::next_collisions(Particle &a, PQ &events)
{
	const vec2d &v = a.v;
	double tx = v.x < 0 ? std::max(0.,  a.p.x-a.r)/-v.x
		  : v.x > 0 ? std::min(w ,w-a.p.x-a.r)/ v.x
		  : INFINITY;
	if (tx < INFINITY)
		events.insert(std::make_shared<C1Ev>(tx, a, true, *this, events));
	double ty = v.y < 0 ? std::max(0.,  a.p.y-a.r)/-v.y
		  : v.y > 0 ? std::min(h ,h-a.p.y-a.r)/ v.y
		  : INFINITY;
	if (ty < INFINITY)
		events.insert(std::make_shared<C1Ev>(ty, a, false, *this, events));

	for (unsigned j=0; j<particles.size(); j++) {
		Particle &b = particles[j];
		auto r = time_to_collision(a, b);
		if (!r.first)
			continue;
		events.insert(std::make_shared<C2Ev>(r.second, a, b, *this, events));
	}
}

std::shared_ptr<Ev> Field::next_collision(Particle &a, PQ &events)
{
	std::shared_ptr<Ev> ev;

	const vec2d &v = a.v;
	double tx = v.x < 0 ? std::max(0.,  a.p.x-a.r)/-v.x
		  : v.x > 0 ? std::min(w ,w-a.p.x-a.r)/ v.x
		  : INFINITY;
	if (tx < INFINITY)
		ev = std::make_shared<C1Ev>(tx, a, true, *this, events);
	double ty = v.y < 0 ? std::max(0.,  a.p.y-a.r)/-v.y
		  : v.y > 0 ? std::min(h ,h-a.p.y-a.r)/ v.y
		  : INFINITY;
	if (ty < INFINITY && (!ev || ty < ev->t))
		ev = std::make_shared<C1Ev>(ty, a, false, *this, events);

	for (unsigned j=0; j<particles.size(); j++) {
		Particle &b = particles[j];
		auto r = time_to_collision(a, b);
		if (!r.first)
			continue;
		if (!ev || r.second < ev->t)
			ev = std::make_shared<C2Ev>(r.second, a, b, *this, events);
	}

	return ev;
}

static constexpr const double EPS = 1e-12;
static double FPS = 60;

#ifdef __unix__
#include <cairo.h>

struct netpbm {
	FILE *out;
	unsigned w, h;
	netpbm(FILE *out) : out(out) {}
	cairo_surface_t * init(unsigned w, unsigned h) noexcept
	{
		this->w = w;
		this->h = h;
		return cairo_image_surface_create(CAIRO_FORMAT_A8, w, h);
	}
	void operator()(cairo_surface_t *sf) const
	{
		cairo_surface_flush(sf);
		unsigned char *data = cairo_image_surface_get_data(sf);
		unsigned s = cairo_image_surface_get_stride(sf);
		fprintf(out, "P5 %u %u 255\n", w, h);
		for (unsigned i=0; i<h; data += s, i++)
			fwrite(data, w, 1, out);
	}
};

struct y4m {
	FILE *out;
	unsigned w, h;
	unsigned char *uv, *line, y_lut[256];
	y4m(FILE *out = stdout, double gamma = 2.2) : out(out)
	{
		for (unsigned i=0; i<256; i++)
			y_lut[i] = 16+std::pow(i/255.0, 1/gamma)*220+.5;
	}
	cairo_surface_t * init(unsigned w, unsigned h)
	{
		this->w = w;
		this->h = h;
		uv = (unsigned char *)malloc(w*h/2);
		memset(uv, 0x80, w*h/2);
		line = (unsigned char *)malloc(w);
		fprintf(out, "YUV4MPEG2 W%u H%u F%u:1 Ip A1:1 C420mpeg2\n", w, h, (unsigned)FPS);
		return cairo_image_surface_create(CAIRO_FORMAT_A8, w, h);
	}
	void operator()(cairo_surface_t *sf) const
	{
		cairo_surface_flush(sf);
		unsigned char *data = cairo_image_surface_get_data(sf);
		unsigned s = cairo_image_surface_get_stride(sf);
		fprintf(out, "FRAME\n");
		for (unsigned i=0; i<h; data += s, i++) {
#if 1
			for (unsigned j=0; j<w; j++)
				line[j] = y_lut[data[j]];
			fwrite(line, w, 1, out);
#else
			fwrite(data, w, 1, out);
#endif
		}
		fwrite(uv, w*h/2, 1, out);
	}
};

#include <cairo-xlib.h>

struct x11 {
	Display *dsp;
	cairo_t *cr;
	cairo_surface_t *xlib_sf;
	x11() : dsp(XOpenDisplay(NULL))
	{
		if (!dsp)
			throw "cannot open display";
	}
	~x11() { if (dsp) XCloseDisplay(dsp); }
	x11(x11 &&o) : dsp(o.dsp) { o.dsp = nullptr; };
	x11 & operator=(const x11 &) = delete;

	cairo_surface_t * init(unsigned w, unsigned h) noexcept
	{
		int screen = DefaultScreen(dsp);
		Drawable da = XCreateSimpleWindow(dsp, DefaultRootWindow(dsp), 0, 0, w, h, 0, 0, 0);
		XSelectInput(dsp, da, ButtonPressMask | KeyPressMask);
		XMapWindow(dsp, da);

		xlib_sf = cairo_xlib_surface_create(dsp, da, DefaultVisual(dsp, screen), w, h);
		cairo_xlib_surface_set_size(xlib_sf, w, h);
		return xlib_sf;
	}

	void operator()(cairo_surface_t *sf) const noexcept
	{
		//cairo_surface_flush(sf);
	}
};

static void draw(cairo_t *cr, const Particle &p) noexcept
{
	cairo_arc(cr, p.p.x, p.p.y, p.r, 0, 2*M_PI);
	cairo_fill(cr);
}

static void draw(cairo_t *cr, const ColoredParticle &p) noexcept
{
	cairo_save(cr);
	cairo_set_source_rgba(cr, p.red, p.green, p.blue, 1.0);
	draw(cr, static_cast<const Particle &>(p));
	cairo_restore(cr);
}

template <typename T>
struct cairo_output {
	const unsigned   w, h;
	T                t;
	cairo_surface_t *sf;
	cairo_t         *cr;

	cairo_output(double c, double _w, double _h, T &&_t = T())
	noexcept
	: w((unsigned)std::ceil(_w*c)), h((unsigned)std::ceil(_h*c)), t(std::move(_t)),
	  sf(t.init(w,h)), cr(cairo_create(sf))
	{
		cairo_scale(cr, c, c);
	}

	~cairo_output() noexcept
	{
		cairo_destroy(cr);
		cairo_surface_destroy(sf);
	}

	template <typename It>
	void draw(It begin, It end) const noexcept
	{
		// cairo_push_group(cr);
#if 1
		cairo_set_source_rgba(cr,0,1,0,0);
		cairo_set_operator(cr, CAIRO_OPERATOR_SOURCE);
		cairo_rectangle(cr, 0, 0, w, h);
		cairo_fill(cr);
#else
		cairo_set_source_rgba(cr,0,0,0,0);
		cairo_set_operator(cr, CAIRO_OPERATOR_CLEAR);
		cairo_rectangle(cr, 0, 0, w, h);
		cairo_paint(cr);
#endif
		cairo_set_operator(cr, CAIRO_OPERATOR_OVER);
		cairo_set_source_rgba(cr,1,1,1,1);
		for (; begin != end; ++begin) {/*
			const Particle &p = *begin;
			cairo_arc(cr, p.p.x, p.p.y, p.r, 0, 2*M_PI);
			cairo_fill(cr);*/
			::draw(cr, *begin);
		}
		// cairo_surface_flush(sf);
		// unsigned char *data = cairo_image_surface_get_data(sf);
		// unsigned s = cairo_image_surface_get_stride(sf);
		//cairo_pop_group_to_source(cr);
		//cairo_paint(cr);
		t(sf);
	}
};
#endif

#include <SDL.h>
#ifdef __unix__
# include <signal.h>
#endif

template <Uint32 which>
struct sdl_subsystem {
	explicit sdl_subsystem(bool sdl_sigint = false)
	{
		SDL_version linked;
		SDL_GetVersion(&linked);
#ifdef __unix__
		struct sigaction int_action, term_action;
		bool hint_worked = false;
		if (!sdl_sigint) {
#if SDL_VERSION_ATLEAST(2,0,4)
			if (linked.major == 2
			    && SDL_VERSIONNUM(linked.major,
			                      linked.minor,
			                      linked.patch
			                     ) >= SDL_VERSIONNUM(2,0,4)
			    && SDL_SetHint(SDL_HINT_NO_SIGNAL_HANDLERS, "1")) {
				hint_worked = true;
			} else
#endif
			{
				sigaction(SIGINT , NULL, &int_action);
				sigaction(SIGTERM, NULL, &term_action);
			}
		}
#endif
		int r = SDL_InitSubSystem(which);
#ifdef __unix__
		if (!sdl_sigint && !hint_worked) {
			sigaction(SIGINT, &int_action , NULL);
			sigaction(SIGINT, &term_action, NULL);
		}
#endif
		if (r < 0)
			throw SDL_GetError();
	}
	~sdl_subsystem() { SDL_QuitSubSystem(which); }
	sdl_subsystem(const sdl_subsystem &) = delete;
	sdl_subsystem & operator=(const sdl_subsystem &) = delete;
};

typedef sdl_subsystem<SDL_INIT_EVERYTHING> sdl;

struct sdl_window {

	sdl_subsystem<SDL_INIT_VIDEO> sdl_video;
	SDL_Window *win;

	sdl_window(
		unsigned w, unsigned h,
		Uint32 window_flags = SDL_WINDOW_ALLOW_HIGHDPI
	) 
	{
		win = SDL_CreateWindow("Bouncing Balls",
		                       SDL_WINDOWPOS_UNDEFINED,
		                       SDL_WINDOWPOS_UNDEFINED,
		                       w, h,
		                       window_flags);
		if (!win)
			throw SDL_GetError();
	}
	sdl_window() noexcept : win(nullptr) {}
	sdl_window(sdl_window &&o) noexcept : win(o.win) { o.win = nullptr; }
	sdl_window & operator=(sdl_window o) noexcept { std::swap(win, o.win); return *this; }

	virtual ~sdl_window() noexcept { if (win) SDL_DestroyWindow(win); }

	virtual void swap() const noexcept { SDL_UpdateWindowSurface(win); }
};


#include <fstream>
#include <string>

#include <GL/glew.h>
#include <GL/gl.h>

struct gl_shader {
	GLuint shader;

	gl_shader(gl_shader &&o) : shader(o.shader) { o.shader = 0; }
	gl_shader & operator=(gl_shader o) { std::swap(shader, o.shader); return *this; }

	explicit gl_shader(GLenum shaderType)
	: shader(glCreateShader(shaderType))
	{
		if (shader == 0)
			throw "error creating GL shader";
	}
	~gl_shader() { if (shader != 0) glDeleteShader(shader); }

	void set_source(const char *path) const
	{
		std::vector<char> src;
		char c;
		for (std::ifstream in(path); (c = in.get()), in;)
			src.push_back(c);
		const GLchar *src_data[] = { src.data(), };
		const GLint   src_lens[] = { (GLint)src.size(), };
		glShaderSource(shader, 1, src_data, src_lens);
	}

	bool compile() const noexcept
	{
		GLint compiled;
		glCompileShader(shader);
		glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
		return compiled == GL_TRUE;
	}

	std::string info_log() const
	{
		std::string r;
		GLint blen;
		GLsizei slen;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &blen);
		if (blen > 0)
			r.resize(blen-1);
		r.reserve(blen);
		glGetShaderInfoLog(shader, blen, &slen, &r[0]);
		r.resize(slen);
		return r;
	}
};

struct gl_program {

	GLuint prog;

	gl_program(gl_program &&o) : prog(o.prog) { o.prog = 0; }
	gl_program & operator=(gl_program o) { std::swap(prog, o.prog); return *this; }

	gl_program() : prog(glCreateProgram())
	{
		if (!prog)
			throw "error creating GL program";
	}
	~gl_program() { if (prog) glDeleteProgram(prog); }

	void attach(const gl_shader &shader) noexcept
	{
		glAttachShader(prog, shader.shader);
	}

	void detach(const gl_shader &shader) noexcept
	{
		glDetachShader(prog, shader.shader);
	}

	bool link() const noexcept
	{
		GLint linked;
		glLinkProgram(prog);
		glGetProgramiv(prog, GL_LINK_STATUS, &linked);
		return linked == GL_TRUE;
	}

	std::string info_log() const
	{
		std::string r;
		GLint blen;
		GLsizei slen;
		glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &blen);
		if (blen > 0)
			r.resize(blen-1);
		r.reserve(blen);
		glGetProgramInfoLog(prog, blen, &slen, &r[0]);
		r.resize(slen);
		return r;
	}
};

#define ARRAY_SIZE(v)	(sizeof(v)/sizeof(*(v)))

static gl_shader gl_shader_prep(GLenum type, const char *src_path)
{
	gl_shader s(type);
	s.set_source(src_path);
	bool compiled = s.compile();
	std::string log = s.info_log();
	const char *stype = "unknown-type";
	switch (type) {
	case GL_VERTEX_SHADER  : stype = "vertex"; break;
	case GL_FRAGMENT_SHADER: stype = "fragment"; break;
	}
	if (log.size() > 0)
		fprintf(stderr, "%s shader compiler log:\n%s", stype, log.c_str());
	if (!compiled) {
		static char buf[128];
		sprintf(buf, "%s shader compilation error", stype);
		throw buf;
	}
	return s;
}

static gl_program gl_program_prep(
	const std::initializer_list<gl_shader> &shaders,
	const std::initializer_list<std::pair<GLuint,const char *>> &attrs
) {
	gl_program p;
	for (const gl_shader &s : shaders)
		p.attach(s);
	for (const auto &a : attrs)
		glBindAttribLocation(p.prog, a.first, a.second);
	bool linked = p.link();
	std::string log = p.info_log();
	if (log.size() > 0)
		fprintf(stderr, "gl program link result:\n%s", log.c_str());
	if (!linked)
		throw "gl program link error";
	return p;
}

struct sdl_gl_window : public sdl_window {

	SDL_GLContext gl_ctx;

	sdl_gl_window() noexcept : gl_ctx(NULL) {}

	sdl_gl_window(
		unsigned w, unsigned h,
		Uint32 window_flags = SDL_WINDOW_ALLOW_HIGHDPI
	)
	: sdl_window(w, h, window_flags | SDL_WINDOW_OPENGL)
	, gl_ctx(SDL_GL_CreateContext(win))
	{
		if (!gl_ctx)
			throw SDL_GetError();
	}

	sdl_gl_window(sdl_gl_window &&o) noexcept
	: sdl_window(std::forward<sdl_window>(o)), gl_ctx(o.gl_ctx)
	{ o.gl_ctx = NULL; }

	~sdl_gl_window() { if (gl_ctx) SDL_GL_DeleteContext(gl_ctx); }

	sdl_gl_window & operator=(sdl_gl_window o) noexcept
	{
		sdl_window::operator=(std::forward<sdl_window>(o));
		std::swap(gl_ctx, o.gl_ctx);
		return *this;
	}

	sdl_gl_window & operator=(sdl_window o) = delete;

	void set_sync_vblank(bool sync_vblank) const noexcept
	{
		if (!sync_vblank)
			SDL_GL_SetSwapInterval(0);
		else if (SDL_GL_SetSwapInterval(-1) == -1)
			SDL_GL_SetSwapInterval(1);
	}

	/* <http://wiki.libsdl.org/SDL_GL_SwapWindow?highlight=(\bCategoryVideo\b)|(CategoryEnum)|(CategoryStruct)>
	 * "On Mac OS X make sure you bind 0 to the draw framebuffer
	 * before swapping the window, otherwise nothing will happen.
	 * See this blog post for more info."
	 * <http://renderingpipeline.com/2012/05/nsopenglcontext-flushbuffer-might-not-do-what-you-think/> */
	void swap() const noexcept override { SDL_GL_SwapWindow(win); }
};

static sdl_gl_window sdl_gl_window_setup(
	unsigned w, unsigned h,
	unsigned gl_major, unsigned gl_minor, SDL_GLprofile profile,
	bool sync_vblank, int msaa, Uint32 window_flags
) {
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, gl_major);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, gl_minor);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, profile);
	if (msaa > 0) {
		SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
		SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, msaa);
	} else {
		SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 0);
		SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 0);
	}
	sdl_gl_window win(w, h, window_flags);
	win.set_sync_vblank(sync_vblank);
	return win;
}

struct sdl_gl_output {

	static constexpr const GLuint DIAM = 0;
	static constexpr const GLuint POS  = 1;
	static constexpr const GLuint COL  = 2;

	sdl_gl_window win;
	GLuint vao, vbo[3];

private:
	void gl_buffer_setup(const Field &f, double c)
	{
		glGenBuffers(ARRAY_SIZE(vbo), vbo);
		size_t n = f.particles.size();

		glBindBuffer(GL_ARRAY_BUFFER, vbo[DIAM]);
		glBufferData(GL_ARRAY_BUFFER, n*sizeof(double), NULL, GL_STATIC_DRAW);
		glVertexAttribPointer(DIAM, 1, GL_DOUBLE, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(DIAM);
		double *d = (double *)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
		for (unsigned i=0; i<n; i++)
			d[i] = f.particles[i].r * c * 2;
		glUnmapBuffer(GL_ARRAY_BUFFER);

		glBindBuffer(GL_ARRAY_BUFFER, vbo[POS]);
		glBufferData(GL_ARRAY_BUFFER, n*2*sizeof(double), NULL, GL_STREAM_DRAW);
		glVertexAttribPointer(POS, 2, GL_DOUBLE, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(POS);

		glBindBuffer(GL_ARRAY_BUFFER, vbo[COL]);
		glBufferData(GL_ARRAY_BUFFER, n*3*sizeof(double), NULL, GL_STREAM_DRAW);
		glVertexAttribPointer(COL, 3, GL_DOUBLE, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(COL);
	}

	void gl_state_setup(
		const Field &f, double c, double _w, double _h,
		const char *vert_shader_path, const char *frag_shader_path
	) {
		glEnable(GL_PROGRAM_POINT_SIZE);
		glEnable(GL_POINT_SPRITE); // needed for r300g to enable gl_PointCoord

		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		gl_buffer_setup(f, c);

		auto shaders = {
			gl_shader_prep(GL_VERTEX_SHADER  , vert_shader_path),
			gl_shader_prep(GL_FRAGMENT_SHADER, frag_shader_path),
		};
		auto attrs = {
			std::make_pair(POS , "in_position"),
			std::make_pair(DIAM, "in_diameter"),
			std::make_pair(COL , "in_color"),
		};
		gl_program p = gl_program_prep(shaders, attrs);
		glUseProgram(p.prog);

		GLfloat sx = (1 - -1)/_w;
		GLfloat sy = (1 - -1)/_h;
		GLfloat transform[4*3] = {
			sx,  0, -1,
			 0, sy, -1,
			 0,  0,  0,
			 0,  0,  1,
		};
		GLint transform_loc = glGetUniformLocation(p.prog, "transform");
		glUniformMatrix3x4fv(transform_loc, 1, GL_TRUE, transform);
	}
public:
	sdl_gl_output(const sdl_gl_output &) = delete;
	sdl_gl_output & operator=(const sdl_gl_output &) = delete;

	sdl_gl_output(
		double c, double _w, double _h,
		Field &f, bool sync_vblank, int msaa,
		Uint32 window_flags = SDL_WINDOW_ALLOW_HIGHDPI,
		const char *vert_shader_path = "bm-vert.glsl",
		const char *frag_shader_path = "bm-frag.glsl"
	)
	: win(sdl_gl_window_setup(std::ceil(c*_w), std::ceil(c*_h),
	                          2, 0, SDL_GL_CONTEXT_PROFILE_CORE,
	                          sync_vblank, msaa, window_flags))
	{
		glewExperimental = GL_TRUE;
		GLenum r = glewInit();
		if (r)
			throw glewGetErrorString(r);

		gl_state_setup(f, c, _w, _h, vert_shader_path, frag_shader_path);
	}

	~sdl_gl_output() noexcept
	{
		glDisableVertexAttribArray(DIAM);
		glDisableVertexAttribArray(POS);
		glDisableVertexAttribArray(COL);
		glDeleteBuffers(ARRAY_SIZE(vbo), vbo);
		glDeleteVertexArrays(1, &vao);
	}

	mutable Uint32 last_draw_at_ticks = 0;
	mutable unsigned nframes = 0;
	mutable Uint32 fps_ticks;

	template <typename It>
	void draw(It begin, It end) const noexcept
	{
		glBindBuffer(GL_ARRAY_BUFFER, vbo[COL]);
		double (*col)[3] = (double (*)[3])glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
		unsigned i = 0;
		for (It it = begin; it != end; ++it) {
			const ColoredParticle &p = *it;
			col[i][0] = p.red;
			col[i][1] = p.green;
			col[i][2] = p.blue;
			i++;
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);

		glBindBuffer(GL_ARRAY_BUFFER, vbo[POS]);
		double (*pos)[2] = (double (*)[2])glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
		i = 0;
		for (It it = begin; it != end; ++it, ++i) {
			const Particle &p = *it;
			pos[i][0] = p.p.x;
			pos[i][1] = p.p.y;
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);

		glClearColor(0,0,0,1);
		glClear(GL_COLOR_BUFFER_BIT);

		glDrawArrays(GL_POINTS, 0, i);

		Uint32 ms_delay = 1000 / FPS;
		Uint32 now = SDL_GetTicks();
		if (now - last_draw_at_ticks + 2 < ms_delay)
			SDL_Delay(ms_delay - (now - last_draw_at_ticks));
		if (now - fps_ticks >= 1000) {
			static char buf[256];
			sprintf(buf, "Bouncing Balls - %.2f FPS", nframes / ((now - fps_ticks) / 1000.0));
			SDL_SetWindowTitle(win.win, buf);
			fps_ticks = now;
			nframes = 0;
		}
		win.swap();
		last_draw_at_ticks = SDL_GetTicks();
		nframes++;
	}
};

template <typename output>
struct FrameEv : public Ev {
	const Field &f;
	PQ &events;
	const output &o;

	FrameEv(double t, const Field &f, PQ &events, const output &o)
	noexcept : Ev(t), f(f), events(events), o(o) {}

	bool handles(const Particle &) const noexcept { return false; }

	void handle() const noexcept
	{
		o.draw(f.particles.cbegin(), f.particles.cend());
		events.insert(std::make_shared<FrameEv>(1.0/FPS, f, events, o));
	}
};

#ifndef M_PI
# define M_PI	acos(-1)
#endif

extern "C"
int main(int argc, char **argv)
{
	int opt;
	unsigned n = 375;
	double scale = 4;
	double aspect = 16.0 / 9.0;
	double w = 250;
	double h = w / aspect;
	double rv = 19;
	double exp_lambda = .75;
	double max_t = INFINITY;
	bool sync_vblank = true;
	int msaa = 0;
	while ((opt = getopt(argc, argv, ":a:cd:e:f:hn:s:St:v:")) != -1)
		switch (opt) {
		case 'a': msaa = atoi(optarg); break;
		case 'c': collision_swap_colors = !collision_swap_colors; break;
		case 'd':
			w = strtof(optarg, &optarg);
			if (*optarg && *optarg != 'x') {
				fprintf(stderr, "illegal parameter for option '-d'\n");
				return 1;
			}
			h = *optarg ? atof(optarg+1) : w / aspect;
			break;
		case 'e': exp_lambda = atof(optarg); break;
		case 'f': FPS = atoi(optarg); break;
		case 'h':
			fprintf(stderr, "usage: %s [-a MSAA] [-c] [-d W[xH]] [-e EXP_LAMBDA] [-f FPS] [-n NUM] [-s SCALE] [-S] [-t MAX_T] [-v MAX_VELOCITY]\n", argv[0]);
			return 0;
		case 'n': n = atoi(optarg); break;
		case 's': scale = atof(optarg); break;
		case 'S': sync_vblank = !sync_vblank; break;
		case 't': max_t = atof(optarg); break;
		case 'v': rv = atof(optarg); break;
		case ':': fprintf(stderr, "option '-%c' needs a parameter\n", optopt); return 1;
		case '?': fprintf(stderr, "unknown option '-%c'\n", optopt); return 1;
		}

	Field f(w, h);

	std::default_random_engine e/*(std::random_device{}())*/;
	using P = std::uniform_real_distribution<>;

	P v(-rv, rv);
	std::exponential_distribution<> r(exp_lambda);

	for (unsigned i=0; i<n; i++) {
		ColoredParticle p;
		p.r = .5+std::min(10.0,r(e)/1);
		p.m = 1.0*(p.r*p.r*M_PI);
		p.v.x = v(e);
		p.v.y = v(e);
		p.set_hsv(P(0,1)(e));
		do {
			p.p.x = P(p.r*2, f.w-p.r*2)(e);
			p.p.y = P(p.r*2, f.h-p.r*2)(e);
		} while (!f.fits(p));
		f.particles.push_back(p);
	}
	try {
	double t = 0;
	//cairo_output<y4m> out(scale, f.w, f.h);
	//cairo_output<x11> out(scale, f.w, f.h);
	sdl_gl_output out(scale, f.w, f.h, f, sync_vblank, msaa);
	PQ events;
	for (Particle &p : f.particles)
		events.insert(f.next_collision(p, events));
	double fps = 1.0 / FPS;
	events.insert(std::make_shared<FrameEv<decltype(out)>>(fps,f,events,out));
	while (t <= max_t) {
		if (events.empty()) {
			std::cerr << "no more collisions\n";
			break;
		}
		double ct = 0;
		for (double tmax = (*events.begin())->t + 0*EPS;
		     !events.empty() && (*events.begin())->t <= tmax;) {
			std::shared_ptr<Ev> ev = *events.begin();
			events.erase(events.begin());
			// std::cerr << "ev: " << events.size() << ", t=" << t+ct << ", " << *ev << std::endl;
			double dt = ev->t - ct;
			f.advance(dt);
			for (const std::shared_ptr<Ev> &e : events)
				e->t -= dt;
			ct = ev->t;
			ev->handle();
			break;
		}
		t += ct;
	}
	} catch (const char *msg) { fprintf(stderr, "FATAL: %s\n", msg); }
}
