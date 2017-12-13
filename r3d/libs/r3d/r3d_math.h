/**
 * r3d -- 3D rendering library
 * author: Andreas Mantler (ands)
 */

#ifndef R3D_MATH_H
#define R3D_MATH_H

#if defined (ARM_MATH_CM4) || (ARM_MATH_CM3) || (ARM_MATH_CM0)
#define R3D_ARM_MATH
#include "arm_math.h"
#else
#include <math.h>
#endif

//#include "core_cm4.h"

#define float_pi		3.1415926535897932f
#define float_pi_over180	0.017453293f

//#pragma anon_unions
typedef struct {
	union {
		float v[2];
		struct {
			float x, y;
		}raw;
		struct {
			float r, g;
		}color;
	}vect;
} vec2_t;
//#pragma anon_unions
typedef struct {
	union {
		float v[3];
		struct {
			float x, y, z;
		}raw;
		struct {
			float r, g, b;
		}color;
		vec2_t xy;
	}vect;
} vec3_t;
//#pragma anon_unions
typedef struct {
	union {
		float v[4];
		struct {
			float x, y, z, w;
		}raw;
		struct {
			float r, g, b, a;
		}color;
		vec2_t xy;
		vec3_t xyz;
		vec3_t rgb;
	}vect;
} vec4_t;
//#pragma anon_unions
typedef struct {
	union {
		float m[16];
		struct {
			float m00, m10, m20, m30,
			      m01, m11, m21, m31,
			      m02, m12, m22, m32,
			      m03, m13, m23, m33;
		}raw;
		vec4_t c[4];
	}mat;
} mat4_t;

#define int_min(a, b) ((a < b) ? a : b)
//static inline int int_min(int a, int b)
//{
//	return ((a < b) ? a : b);
//}
#define int_max(a, b) ((a < b) ? b : a)
//static inline int int_max(int a, int b)
//{
//	return ((a < b) ? b : a);
//}
#define int_clamp(value, min, max) (value < min ? min : (value < max ? value : max))
//static inline int int_clamp(int value, int min, int max)
//{
//	return value < min ? min : (value < max ? value : max);
//}
#define float_min(a, b) int_min(a, b)
//static inline float float_min(float a, float b)
//{
//	return a < b ? a : b;
//}
#define float_max(a, b) int_max(a, b)
//static inline float float_max(float a, float b)
//{
//	return a < b ? b : a;
//}
#define float_clamp(value, min, max) int_clamp(value, min, max)
//static inline float float_clamp(float value, float min, float max)
//{
//	return value < min ? min : (value < max ? value : max);
//}

static __INLINE vec2_t vec2(float x, float y)
{
	vec2_t r;
	r.vect.raw.x = x;
	r.vect.raw.y = y;
	return r;
}
static __INLINE vec2_t vec2_add(vec2_t v0, vec2_t v1)
{
	vec2_t r;
	r.vect.raw.x = v0.vect.raw.x + v1.vect.raw.x;
	r.vect.raw.y = v0.vect.raw.y + v1.vect.raw.y;
	return r;
}
static __INLINE vec2_t vec2_sub(vec2_t v0, vec2_t v1)
{
	vec2_t r;
	r.vect.raw.x = v0.vect.raw.x - v1.vect.raw.x;
	r.vect.raw.y = v0.vect.raw.y - v1.vect.raw.y;
	return r;
}
static __INLINE vec2_t vec2_mul(vec2_t v, float f)
{
	vec2_t r;
	r.vect.raw.x = v.vect.raw.x * f;
	r.vect.raw.y = v.vect.raw.y * f;
	return r;
}
static __INLINE vec2_t vec2_div(vec2_t v, float f)
{
	float fi = 1.0f / f;
	vec2_t r;
	r.vect.raw.x = v.vect.raw.x * fi;
	r.vect.raw.y = v.vect.raw.y * fi;
	return r;
}
static __INLINE float vec2_dot(vec2_t v0, vec2_t v1)
{
	return v0.vect.raw.x * v1.vect.raw.x + v0.vect.raw.y * v1.vect.raw.y;
}
static __INLINE float vec2_cross(vec2_t v0, vec2_t v1)
{
	return v0.vect.raw.x * v1.vect.raw.y - v0.vect.raw.y * v1.vect.raw.x;
}
static __INLINE float vec2_length(vec2_t v)
{
	return sqrtf(v.vect.raw.x * v.vect.raw.x + v.vect.raw.y * v.vect.raw.y);
}
static __INLINE vec2_t vec2_normalize(vec2_t v)
{
	float fi = 1.0f / sqrtf(v.vect.raw.x * v.vect.raw.x + v.vect.raw.y * v.vect.raw.y);
	vec2_t r;
	r.vect.raw.x = v.vect.raw.x * fi;
	r.vect.raw.y = v.vect.raw.y * fi;
	return r;
}

static __INLINE vec3_t vec3(float x, float y, float z)
{
	vec3_t r;
	r.vect.raw.x = x;
	r.vect.raw.y = y;
	r.vect.raw.z = z;
	return r;
}
static __INLINE vec3_t vec3_2(vec2_t xy, float z)
{
	vec3_t r;
	r.vect.raw.x = xy.vect.raw.x;
	r.vect.raw.y = xy.vect.raw.y;
	r.vect.raw.z = z;
	return r;
}
static __INLINE vec3_t vec3_add(vec3_t v0, vec3_t v1)
{
	vec3_t r;
	r.vect.raw.x = v0.vect.raw.x + v1.vect.raw.x;
	r.vect.raw.y = v0.vect.raw.y + v1.vect.raw.y;
	r.vect.raw.z = v0.vect.raw.z + v1.vect.raw.z;
	return r;
}
static __INLINE vec3_t vec3_sub(vec3_t v0, vec3_t v1)
{
	vec3_t r;
	r.vect.raw.x = v0.vect.raw.x - v1.vect.raw.x;
	r.vect.raw.y = v0.vect.raw.y - v1.vect.raw.y;
	r.vect.raw.z = v0.vect.raw.z - v1.vect.raw.z;
	return r;
}
static __INLINE vec3_t vec3_mul(vec3_t v, float f)
{
	vec3_t r;
	r.vect.raw.x = v.vect.raw.x * f;
	r.vect.raw.y = v.vect.raw.y * f;
	r.vect.raw.z = v.vect.raw.z * f;
	return r;
}
static __INLINE vec3_t vec3_div(vec3_t v, float f)
{
	float fi = 1.0f / f;
	vec3_t r;
	r.vect.raw.x = v.vect.raw.x * fi;
	r.vect.raw.y = v.vect.raw.y * fi;
	r.vect.raw.z = v.vect.raw.z * fi;
	return r;
}
static __INLINE float vec3_dot(vec3_t v0, vec3_t v1)
{
	return v0.vect.raw.x * v1.vect.raw.x + v0.vect.raw.y * v1.vect.raw.y + v0.vect.raw.z * v1.vect.raw.z;
}
static __INLINE vec3_t vec3_cross(vec3_t v0, vec3_t v1)
{
	vec3_t r;
	r.vect.raw.x = v0.vect.raw.y * v1.vect.raw.z - v0.vect.raw.z * v1.vect.raw.y;
	r.vect.raw.y = v0.vect.raw.z * v1.vect.raw.x - v0.vect.raw.x * v1.vect.raw.z;
	r.vect.raw.z = v0.vect.raw.x * v1.vect.raw.y - v0.vect.raw.y * v1.vect.raw.x;
	return r;
}
static __INLINE float vec3_length(vec3_t v)
{
	return sqrtf(v.vect.raw.x * v.vect.raw.x + v.vect.raw.y * v.vect.raw.y + v.vect.raw.z * v.vect.raw.z);
}
static __INLINE vec3_t vec3_normalize(vec3_t v)
{
	float fi = 1.0f / sqrtf(v.vect.raw.x * v.vect.raw.x + v.vect.raw.y * v.vect.raw.y + v.vect.raw.z * v.vect.raw.z);
	vec3_t r;
	r.vect.raw.x = v.vect.raw.x * fi;
	r.vect.raw.y = v.vect.raw.y * fi;
	r.vect.raw.z = v.vect.raw.z * fi;
	return r;
}

static __INLINE vec4_t vec4(float x, float y, float z, float w)
{
	vec4_t r;
	r.vect.raw.x = x;
	r.vect.raw.y = y;
	r.vect.raw.z = z;
	r.vect.raw.w = w;
	return r;
}
static __INLINE vec4_t vec4_2(const vec2_t xy, const float z, const float w)
{
	vec4_t r;
	r.vect.raw.x = xy.vect.raw.x;
	r.vect.raw.y = xy.vect.raw.y;
	r.vect.raw.z = z;
	r.vect.raw.w = w;
	return r;
}
static __INLINE vec4_t vec4_3(vec3_t xyz, float w)
{
	vec4_t r;
	r.vect.raw.x = xyz.vect.raw.x;
	r.vect.raw.y = xyz.vect.raw.y;
	r.vect.raw.z = xyz.vect.raw.z;
	r.vect.raw.w = w;
	return r;
}
static __INLINE vec4_t vec4_add(vec4_t v0, vec4_t v1)
{
	vec4_t r;
	r.vect.raw.x = v0.vect.raw.x + v1.vect.raw.x;
	r.vect.raw.y = v0.vect.raw.y + v1.vect.raw.y;
	r.vect.raw.z = v0.vect.raw.z + v1.vect.raw.z;
	r.vect.raw.w = v0.vect.raw.w + v1.vect.raw.w;
	return r;
}
static __INLINE vec4_t vec4_sub(vec4_t v0, vec4_t v1)
{
	vec4_t r;
	r.vect.raw.x = v0.vect.raw.x - v1.vect.raw.x;
	r.vect.raw.y = v0.vect.raw.y - v1.vect.raw.y;
	r.vect.raw.z = v0.vect.raw.z - v1.vect.raw.z;
	r.vect.raw.w = v0.vect.raw.w - v1.vect.raw.w;
	return r;
}
static __INLINE vec4_t vec4_mul(vec4_t v, float f)
{
	vec4_t r;
	r.vect.raw.x = v.vect.raw.x * f;
	r.vect.raw.y = v.vect.raw.y * f;
	r.vect.raw.z = v.vect.raw.z * f;
	r.vect.raw.w = v.vect.raw.w * f;
	return r;
}
static __INLINE vec4_t vec4_div(vec4_t v, float f)
{
	float fi = 1.0f / f;
	vec4_t r;
	r.vect.raw.x = v.vect.raw.x * fi;
	r.vect.raw.y = v.vect.raw.y * fi;
	r.vect.raw.z = v.vect.raw.z * fi;
	r.vect.raw.w = v.vect.raw.w * fi;
	return r;
}
static __INLINE float vec4_dot(vec4_t v0, vec4_t v1)
{
	return v0.vect.raw.x * v1.vect.raw.x + v0.vect.raw.y * v1.vect.raw.y + v0.vect.raw.z * v1.vect.raw.z + v0.vect.raw.w * v1.vect.raw.w;
}
static __INLINE vec4_t vec4_cross(vec4_t v0, vec4_t v1)
{
	return vec4_3(vec3_cross(v0.vect.xyz, v1.vect.xyz), (v0.vect.raw.w + v1.vect.raw.w) * 0.5f);    // USUALLY FAIL
}
static __INLINE float vec4_length(vec4_t v)
{
	return sqrtf(v.vect.raw.x * v.vect.raw.x + v.vect.raw.y * v.vect.raw.y + v.vect.raw.z * v.vect.raw.z + v.vect.raw.w * v.vect.raw.w);
}
static __INLINE vec4_t vec4_normalize(vec4_t v)
{
	float fi = 1.0f / sqrtf(v.vect.raw.x * v.vect.raw.x + v.vect.raw.y * v.vect.raw.y + v.vect.raw.z * v.vect.raw.z + v.vect.raw.w * v.vect.raw.w);
	vec4_t r;
	r.vect.raw.x = v.vect.raw.x * fi;
	r.vect.raw.y = v.vect.raw.y * fi;
	r.vect.raw.z = v.vect.raw.z * fi;
	r.vect.raw.w = v.vect.raw.w * fi;
	return r;
}

static __INLINE mat4_t mat4(float m00, float m10, float m20, float m30, float m01, float m11, float m21, float m31, float m02, float m12, float m22, float m32, float m03, float m13, float m23, float m33)
{
	mat4_t m;
	m.mat.raw.m00 = m00; m.mat.raw.m10 = m10; m.mat.raw.m20 = m20; m.mat.raw.m30 = m30;
	m.mat.raw.m01 = m01; m.mat.raw.m11 = m11; m.mat.raw.m21 = m21; m.mat.raw.m31 = m31;
	m.mat.raw.m02 = m02; m.mat.raw.m12 = m12; m.mat.raw.m22 = m22; m.mat.raw.m32 = m32;
	m.mat.raw.m03 = m03; m.mat.raw.m13 = m13; m.mat.raw.m23 = m23; m.mat.raw.m33 = m33;
	return m;
}

static __INLINE mat4_t mat4_identity()
{
	mat4_t m;
	m.mat.raw.m00 = 1.0f; m.mat.raw.m10 = 0.0f; m.mat.raw.m20 = 0.0f; m.mat.raw.m30 = 0.0f;
	m.mat.raw.m01 = 0.0f; m.mat.raw.m11 = 1.0f; m.mat.raw.m21 = 0.0f; m.mat.raw.m31 = 0.0f;
	m.mat.raw.m02 = 0.0f; m.mat.raw.m12 = 0.0f; m.mat.raw.m22 = 1.0f; m.mat.raw.m32 = 0.0f;
	m.mat.raw.m03 = 0.0f; m.mat.raw.m13 = 0.0f; m.mat.raw.m23 = 0.0f; m.mat.raw.m33 = 1.0f;
	return m;
}

static __INLINE mat4_t mat4_ortho(float left, float right,
                                float bottom, float top,
                                float near, float far)
{
	float rl = right - left;
	float tx = - (right + left) / rl;
	float tb = top - bottom;
	float ty = - (top + bottom) / tb;
	float fn = far - near;
	float tz = - (far + near) / fn;
	mat4_t m;
	m.mat.raw.m00 = 2.0f / rl; m.mat.raw.m10 = 0.0f; m.mat.raw.m20 = 0.0f; m.mat.raw.m30 = 0.0f;
	m.mat.raw.m01 = 0.0f; m.mat.raw.m11 = 2.0f / tb; m.mat.raw.m21 = 0.0f; m.mat.raw.m31 = 0.0f;
	m.mat.raw.m02 = 0.0f; m.mat.raw.m12 = 0.0f; m.mat.raw.m22 = -2.0f / fn; m.mat.raw.m32 = 0.0f;
	m.mat.raw.m03 = tx; m.mat.raw.m13 = ty; m.mat.raw.m23 = tz; m.mat.raw.m33 = 1.0f;
	return m;
}

static __INLINE mat4_t mat4_perspective(float fovy, float aspect,
                                      float near, float far)
{
	float a = fovy * float_pi_over180 * 0.5f;
#ifdef R3D_ARM_MATH
	float f = arm_cos_f32(a) / arm_sin_f32(a);
#else
	float f = 1.0f / tanf(a);
#endif
	float nfi = 1.0f / (near - far);
	mat4_t m;
	m.mat.raw.m00 = f / aspect; m.mat.raw.m10 = 0.0f; m.mat.raw.m20 = 0.0f; m.mat.raw.m30 = 0.0f;
	m.mat.raw.m01 = 0.0f; m.mat.raw.m11 = f; m.mat.raw.m21 = 0.0f; m.mat.raw.m31 = 0.0f;
	m.mat.raw.m02 = 0.0f; m.mat.raw.m12 = 0.0f; m.mat.raw.m22 = (far + near) * nfi; m.mat.raw.m32 = -1.0f;
	m.mat.raw.m03 = 0.0f; m.mat.raw.m13 = 0.0f; m.mat.raw.m23 = 2.0f * far * near * nfi; m.mat.raw.m33 = 0.0f;
	return m;
}

static __INLINE mat4_t mat4_mul(mat4_t m0, mat4_t m1)
{
	mat4_t m;
	m.mat.raw.m00 = m0.mat.raw.m00*m1.mat.raw.m00 + m0.mat.raw.m01*m1.mat.raw.m10 + m0.mat.raw.m02*m1.mat.raw.m20 + m0.mat.raw.m03 * m1.mat.raw.m30;
	m.mat.raw.m10 = m0.mat.raw.m10*m1.mat.raw.m00 + m0.mat.raw.m11*m1.mat.raw.m10 + m0.mat.raw.m12*m1.mat.raw.m20 + m0.mat.raw.m13 * m1.mat.raw.m30;
	m.mat.raw.m20 = m0.mat.raw.m20*m1.mat.raw.m00 + m0.mat.raw.m21*m1.mat.raw.m10 + m0.mat.raw.m22*m1.mat.raw.m20 + m0.mat.raw.m23 * m1.mat.raw.m30;
	m.mat.raw.m30 = m0.mat.raw.m30*m1.mat.raw.m00 + m0.mat.raw.m31*m1.mat.raw.m10 + m0.mat.raw.m32*m1.mat.raw.m20 + m0.mat.raw.m33 * m1.mat.raw.m30;
                  
	m.mat.raw.m01 = m0.mat.raw.m00*m1.mat.raw.m01 + m0.mat.raw.m01*m1.mat.raw.m11 + m0.mat.raw.m02*m1.mat.raw.m21 + m0.mat.raw.m03 * m1.mat.raw.m31;
	m.mat.raw.m11 = m0.mat.raw.m10*m1.mat.raw.m01 + m0.mat.raw.m11*m1.mat.raw.m11 + m0.mat.raw.m12*m1.mat.raw.m21 + m0.mat.raw.m13 * m1.mat.raw.m31;
	m.mat.raw.m21 = m0.mat.raw.m20*m1.mat.raw.m01 + m0.mat.raw.m21*m1.mat.raw.m11 + m0.mat.raw.m22*m1.mat.raw.m21 + m0.mat.raw.m23 * m1.mat.raw.m31;
	m.mat.raw.m31 = m0.mat.raw.m30*m1.mat.raw.m01 + m0.mat.raw.m31*m1.mat.raw.m11 + m0.mat.raw.m32*m1.mat.raw.m21 + m0.mat.raw.m33 * m1.mat.raw.m31;
                  
	m.mat.raw.m02 = m0.mat.raw.m00*m1.mat.raw.m02 + m0.mat.raw.m01*m1.mat.raw.m12 + m0.mat.raw.m02*m1.mat.raw.m22 + m0.mat.raw.m03 * m1.mat.raw.m32;
	m.mat.raw.m12 = m0.mat.raw.m10*m1.mat.raw.m02 + m0.mat.raw.m11*m1.mat.raw.m12 + m0.mat.raw.m12*m1.mat.raw.m22 + m0.mat.raw.m13 * m1.mat.raw.m32;
	m.mat.raw.m22 = m0.mat.raw.m20*m1.mat.raw.m02 + m0.mat.raw.m21*m1.mat.raw.m12 + m0.mat.raw.m22*m1.mat.raw.m22 + m0.mat.raw.m23 * m1.mat.raw.m32;
	m.mat.raw.m32 = m0.mat.raw.m30*m1.mat.raw.m02 + m0.mat.raw.m31*m1.mat.raw.m12 + m0.mat.raw.m32*m1.mat.raw.m22 + m0.mat.raw.m33 * m1.mat.raw.m32;
                  
	m.mat.raw.m03 = m0.mat.raw.m00*m1.mat.raw.m03 + m0.mat.raw.m01*m1.mat.raw.m13 + m0.mat.raw.m02*m1.mat.raw.m23 + m0.mat.raw.m03 * m1.mat.raw.m33;
	m.mat.raw.m13 = m0.mat.raw.m10*m1.mat.raw.m03 + m0.mat.raw.m11*m1.mat.raw.m13 + m0.mat.raw.m12*m1.mat.raw.m23 + m0.mat.raw.m13 * m1.mat.raw.m33;
	m.mat.raw.m23 = m0.mat.raw.m20*m1.mat.raw.m03 + m0.mat.raw.m21*m1.mat.raw.m13 + m0.mat.raw.m22*m1.mat.raw.m23 + m0.mat.raw.m23 * m1.mat.raw.m33;
	m.mat.raw.m33 = m0.mat.raw.m30*m1.mat.raw.m03 + m0.mat.raw.m31*m1.mat.raw.m13 + m0.mat.raw.m32*m1.mat.raw.m23 + m0.mat.raw.m33 * m1.mat.raw.m33;
	return m;
}

static __INLINE mat4_t mat4_translation(vec3_t v)
{
	mat4_t m;
	m.mat.raw.m00 = 1.0f; m.mat.raw.m10 = 0.0f; m.mat.raw.m20 = 0.0f; m.mat.raw.m30 = 0.0f;
	m.mat.raw.m01 = 0.0f; m.mat.raw.m11 = 1.0f; m.mat.raw.m21 = 0.0f; m.mat.raw.m31 = 0.0f;
	m.mat.raw.m02 = 0.0f; m.mat.raw.m12 = 0.0f; m.mat.raw.m22 = 1.0f; m.mat.raw.m32 = 0.0f;
	m.mat.raw.m03 = v.vect.raw.x; m.mat.raw.m13 = v.vect.raw.y; m.mat.raw.m23 = v.vect.raw.z; m.mat.raw.m33 = 1.0f;
	return m;
}

static __INLINE mat4_t mat4_lookat(vec3_t eye, vec3_t center, vec3_t up)
{
	vec3_t f = vec3_normalize(vec3_sub(center, eye));
	vec3_t s = vec3_cross(f, vec3_normalize(up));
	vec3_t u = vec3_cross(vec3_normalize(s), f);

	mat4_t m;
	m.mat.raw.m00 = s.vect.raw.x; m.mat.raw.m10 = u.vect.raw.x; m.mat.raw.m20 = -f.vect.raw.x; m.mat.raw.m30 = 0.0f;
	m.mat.raw.m01 = s.vect.raw.y; m.mat.raw.m11 = u.vect.raw.y; m.mat.raw.m21 = -f.vect.raw.y; m.mat.raw.m31 = 0.0f;
	m.mat.raw.m02 = s.vect.raw.z; m.mat.raw.m12 = u.vect.raw.z; m.mat.raw.m22 = -f.vect.raw.z; m.mat.raw.m32 = 0.0f;
	m.mat.raw.m03 = 0.0f; m.mat.raw.m13 = 0.0f; m.mat.raw.m23 = 0.0f; m.mat.raw.m33 = 1.0f;
	return mat4_mul(m, mat4_translation(vec3_mul(eye, -1.0f)));
}

static __INLINE mat4_t mat4_scaling(vec3_t v)
{
	mat4_t m;
	m.mat.raw.m00 = v.vect.raw.x; m.mat.raw.m10 = 0.0f; m.mat.raw.m20 = 0.0f; m.mat.raw.m30 = 0.0f;
	m.mat.raw.m01 = 0.0f; m.mat.raw.m11 = v.vect.raw.y; m.mat.raw.m21 = 0.0f; m.mat.raw.m31 = 0.0f;
	m.mat.raw.m02 = 0.0f; m.mat.raw.m12 = 0.0f; m.mat.raw.m22 = v.vect.raw.z; m.mat.raw.m32 = 0.0f;
	m.mat.raw.m03 = 0.0f; m.mat.raw.m13 = 0.0f; m.mat.raw.m23 = 0.0f; m.mat.raw.m33 = 1.0f;
	return m;
}

static __INLINE mat4_t mat4_rotation(float angle, vec3_t axis)
{
#ifdef R3D_ARM_MATH
	float c, s;
	const float i360 = 1.0f / 360.0f;
	int a = (int)(angle * i360);
	if (angle < 0.0f) a--;
	arm_sin_cos_f32(angle - (a * 360.0f) - 180.0f, &c, &s); // accepts [-180, 180]
	// = arm_sin_cos_f32(fmod(angle, 360) - 180, &c, &s)
	c *= -1.0f;
	s *= -1.0f; // we're 180° off. correct this
#else
	angle *= float_pi_over180;
	float c = cosf(angle);
	float s = sinf(angle);
#endif
	float c2 = 1.0f - c;
	axis = vec3_normalize(axis);
	float x = axis.vect.raw.x;
	float y = axis.vect.raw.y;
	float z = axis.vect.raw.z;

	mat4_t m;
	m.mat.raw.m00 = x*x*c2 + c; m.mat.raw.m10 = y*x*c2 + z * s; m.mat.raw.m20 = x*z*c2 - y * s; m.mat.raw.m30 = 0.0f;
	m.mat.raw.m01 = x*y*c2 - z * s; m.mat.raw.m11 = y*y*c2 + c; m.mat.raw.m21 = y*z*c2 + x * s; m.mat.raw.m31 = 0.0f;
	m.mat.raw.m02 = x*z*c2 + y * s; m.mat.raw.m12 = y*z*c2 - x * s; m.mat.raw.m22 = z*z*c2 + c; m.mat.raw.m32 = 0.0f;
	m.mat.raw.m03 = 0.0f; m.mat.raw.m13 = 0.0f; m.mat.raw.m23 = 0.0f; m.mat.raw.m33 = 1.0f;
	return m;
}

static __INLINE vec4_t mat4_transform(mat4_t m, vec4_t v)
{
	vec4_t r;
	r.vect.raw.x = m.mat.raw.m00*v.vect.raw.x + m.mat.raw.m01*v.vect.raw.y + m.mat.raw.m02*v.vect.raw.z + m.mat.raw.m03 * v.vect.raw.w;
	r.vect.raw.y = m.mat.raw.m10*v.vect.raw.x + m.mat.raw.m11*v.vect.raw.y + m.mat.raw.m12*v.vect.raw.z + m.mat.raw.m13 * v.vect.raw.w;
	r.vect.raw.z = m.mat.raw.m20*v.vect.raw.x + m.mat.raw.m21*v.vect.raw.y + m.mat.raw.m22*v.vect.raw.z + m.mat.raw.m23 * v.vect.raw.w;
	r.vect.raw.w = m.mat.raw.m30*v.vect.raw.x + m.mat.raw.m31*v.vect.raw.y + m.mat.raw.m32*v.vect.raw.z + m.mat.raw.m33 * v.vect.raw.w;
	return r;
}

static __INLINE vec3_t mat4_transform_position(mat4_t m, vec3_t v)
{
	float fi = 1.0f / (m.mat.raw.m30 * v.vect.raw.x + m.mat.raw.m31 * v.vect.raw.y + m.mat.raw.m32 * v.vect.raw.z + m.mat.raw.m33);
	vec3_t r;
	r.vect.raw.x = (m.mat.raw.m00 * v.vect.raw.x + m.mat.raw.m01 * v.vect.raw.y + m.mat.raw.m02 * v.vect.raw.z + m.mat.raw.m03) * fi;
	r.vect.raw.y = (m.mat.raw.m10 * v.vect.raw.x + m.mat.raw.m11 * v.vect.raw.y + m.mat.raw.m12 * v.vect.raw.z + m.mat.raw.m13) * fi;
	r.vect.raw.z = (m.mat.raw.m20 * v.vect.raw.x + m.mat.raw.m21 * v.vect.raw.y + m.mat.raw.m22 * v.vect.raw.z + m.mat.raw.m23) * fi;
	return r;
}

static __INLINE vec3_t mat4_transform_vector(mat4_t m, vec3_t v)
{
	vec3_t r;
	r.vect.raw.x = m.mat.raw.m00*v.vect.raw.x + m.mat.raw.m01*v.vect.raw.y + m.mat.raw.m02 * v.vect.raw.z;
	r.vect.raw.y = m.mat.raw.m10*v.vect.raw.x + m.mat.raw.m11*v.vect.raw.y + m.mat.raw.m12 * v.vect.raw.z;
	r.vect.raw.z = m.mat.raw.m20*v.vect.raw.x + m.mat.raw.m21*v.vect.raw.y + m.mat.raw.m22 * v.vect.raw.z;
	return r;
}

static __INLINE mat4_t mat4_invert(mat4_t m)
{
	mat4_t mi;
	mi.mat.raw.m00 = m.mat.raw.m11*m.mat.raw.m22*m.mat.raw.m33 + m.mat.raw.m12*m.mat.raw.m23*m.mat.raw.m31 + m.mat.raw.m13*m.mat.raw.m21*m.mat.raw.m32 - m.mat.raw.m11*m.mat.raw.m23*m.mat.raw.m32 - m.mat.raw.m12*m.mat.raw.m21*m.mat.raw.m33 - m.mat.raw.m13*m.mat.raw.m22 * m.mat.raw.m31;
	mi.mat.raw.m10 = m.mat.raw.m10*m.mat.raw.m23*m.mat.raw.m32 + m.mat.raw.m12*m.mat.raw.m20*m.mat.raw.m33 + m.mat.raw.m13*m.mat.raw.m22*m.mat.raw.m30 - m.mat.raw.m10*m.mat.raw.m22*m.mat.raw.m33 - m.mat.raw.m12*m.mat.raw.m23*m.mat.raw.m30 - m.mat.raw.m13*m.mat.raw.m20 * m.mat.raw.m32;
	mi.mat.raw.m20 = m.mat.raw.m10*m.mat.raw.m21*m.mat.raw.m33 + m.mat.raw.m11*m.mat.raw.m23*m.mat.raw.m30 + m.mat.raw.m13*m.mat.raw.m20*m.mat.raw.m31 - m.mat.raw.m10*m.mat.raw.m23*m.mat.raw.m31 - m.mat.raw.m11*m.mat.raw.m20*m.mat.raw.m33 - m.mat.raw.m13*m.mat.raw.m21 * m.mat.raw.m30;
	mi.mat.raw.m30 = m.mat.raw.m10*m.mat.raw.m22*m.mat.raw.m31 + m.mat.raw.m11*m.mat.raw.m20*m.mat.raw.m32 + m.mat.raw.m12*m.mat.raw.m21*m.mat.raw.m30 - m.mat.raw.m10*m.mat.raw.m21*m.mat.raw.m32 - m.mat.raw.m11*m.mat.raw.m22*m.mat.raw.m30 - m.mat.raw.m12*m.mat.raw.m20 * m.mat.raw.m31;

	mi.mat.raw.m01 = m.mat.raw.m01*m.mat.raw.m23*m.mat.raw.m32 + m.mat.raw.m02*m.mat.raw.m21*m.mat.raw.m33 + m.mat.raw.m03*m.mat.raw.m22*m.mat.raw.m31 - m.mat.raw.m01*m.mat.raw.m22*m.mat.raw.m33 - m.mat.raw.m02*m.mat.raw.m23*m.mat.raw.m31 - m.mat.raw.m03*m.mat.raw.m21 * m.mat.raw.m32;
	mi.mat.raw.m11 = m.mat.raw.m00*m.mat.raw.m22*m.mat.raw.m33 + m.mat.raw.m02*m.mat.raw.m23*m.mat.raw.m30 + m.mat.raw.m03*m.mat.raw.m20*m.mat.raw.m32 - m.mat.raw.m00*m.mat.raw.m23*m.mat.raw.m32 - m.mat.raw.m02*m.mat.raw.m20*m.mat.raw.m33 - m.mat.raw.m03*m.mat.raw.m22 * m.mat.raw.m30;
	mi.mat.raw.m21 = m.mat.raw.m00*m.mat.raw.m23*m.mat.raw.m31 + m.mat.raw.m01*m.mat.raw.m20*m.mat.raw.m33 + m.mat.raw.m03*m.mat.raw.m31*m.mat.raw.m30 - m.mat.raw.m00*m.mat.raw.m21*m.mat.raw.m33 - m.mat.raw.m01*m.mat.raw.m23*m.mat.raw.m30 - m.mat.raw.m03*m.mat.raw.m20 * m.mat.raw.m31;
	mi.mat.raw.m31 = m.mat.raw.m00*m.mat.raw.m21*m.mat.raw.m32 + m.mat.raw.m01*m.mat.raw.m22*m.mat.raw.m30 + m.mat.raw.m02*m.mat.raw.m20*m.mat.raw.m31 - m.mat.raw.m00*m.mat.raw.m22*m.mat.raw.m31 - m.mat.raw.m01*m.mat.raw.m20*m.mat.raw.m32 - m.mat.raw.m02*m.mat.raw.m21 * m.mat.raw.m30;

	mi.mat.raw.m02 = m.mat.raw.m01*m.mat.raw.m12*m.mat.raw.m33 + m.mat.raw.m02*m.mat.raw.m13*m.mat.raw.m31 + m.mat.raw.m03*m.mat.raw.m11*m.mat.raw.m32 - m.mat.raw.m01*m.mat.raw.m13*m.mat.raw.m32 - m.mat.raw.m02*m.mat.raw.m11*m.mat.raw.m33 - m.mat.raw.m03*m.mat.raw.m12 * m.mat.raw.m31;
	mi.mat.raw.m12 = m.mat.raw.m00*m.mat.raw.m13*m.mat.raw.m32 + m.mat.raw.m02*m.mat.raw.m10*m.mat.raw.m33 + m.mat.raw.m03*m.mat.raw.m12*m.mat.raw.m30 - m.mat.raw.m00*m.mat.raw.m12*m.mat.raw.m33 - m.mat.raw.m02*m.mat.raw.m13*m.mat.raw.m30 - m.mat.raw.m03*m.mat.raw.m10 * m.mat.raw.m32;
	mi.mat.raw.m22 = m.mat.raw.m00*m.mat.raw.m11*m.mat.raw.m33 + m.mat.raw.m01*m.mat.raw.m13*m.mat.raw.m30 + m.mat.raw.m03*m.mat.raw.m10*m.mat.raw.m31 - m.mat.raw.m00*m.mat.raw.m13*m.mat.raw.m31 - m.mat.raw.m01*m.mat.raw.m10*m.mat.raw.m33 - m.mat.raw.m03*m.mat.raw.m11 * m.mat.raw.m30;
	mi.mat.raw.m32 = m.mat.raw.m00*m.mat.raw.m12*m.mat.raw.m31 + m.mat.raw.m01*m.mat.raw.m10*m.mat.raw.m32 + m.mat.raw.m02*m.mat.raw.m11*m.mat.raw.m30 - m.mat.raw.m00*m.mat.raw.m11*m.mat.raw.m32 - m.mat.raw.m01*m.mat.raw.m12*m.mat.raw.m30 - m.mat.raw.m02*m.mat.raw.m10 * m.mat.raw.m31;

	mi.mat.raw.m03 = m.mat.raw.m01*m.mat.raw.m13*m.mat.raw.m22 + m.mat.raw.m02*m.mat.raw.m11*m.mat.raw.m23 + m.mat.raw.m03*m.mat.raw.m12*m.mat.raw.m21 - m.mat.raw.m01*m.mat.raw.m12*m.mat.raw.m23 - m.mat.raw.m02*m.mat.raw.m13*m.mat.raw.m21 - m.mat.raw.m03*m.mat.raw.m11 * m.mat.raw.m22;
	mi.mat.raw.m13 = m.mat.raw.m00*m.mat.raw.m12*m.mat.raw.m23 + m.mat.raw.m02*m.mat.raw.m13*m.mat.raw.m20 + m.mat.raw.m03*m.mat.raw.m10*m.mat.raw.m22 - m.mat.raw.m00*m.mat.raw.m13*m.mat.raw.m22 - m.mat.raw.m02*m.mat.raw.m10*m.mat.raw.m23 - m.mat.raw.m03*m.mat.raw.m12 * m.mat.raw.m20;
	mi.mat.raw.m23 = m.mat.raw.m00*m.mat.raw.m13*m.mat.raw.m21 + m.mat.raw.m01*m.mat.raw.m10*m.mat.raw.m23 + m.mat.raw.m03*m.mat.raw.m11*m.mat.raw.m20 - m.mat.raw.m00*m.mat.raw.m11*m.mat.raw.m23 - m.mat.raw.m01*m.mat.raw.m13*m.mat.raw.m20 - m.mat.raw.m03*m.mat.raw.m10 * m.mat.raw.m21;
	mi.mat.raw.m33 = m.mat.raw.m00*m.mat.raw.m11*m.mat.raw.m22 + m.mat.raw.m01*m.mat.raw.m12*m.mat.raw.m20 + m.mat.raw.m02*m.mat.raw.m10*m.mat.raw.m21 - m.mat.raw.m00*m.mat.raw.m12*m.mat.raw.m21 - m.mat.raw.m01*m.mat.raw.m10*m.mat.raw.m22 - m.mat.raw.m02*m.mat.raw.m11 * m.mat.raw.m20;
	return mi;
}

static __INLINE mat4_t mat4_transpose(mat4_t m)
{
	mat4_t mt;
	mt.mat.raw.m00 = m.mat.raw.m00; mt.mat.raw.m10 = m.mat.raw.m01; mt.mat.raw.m20 = m.mat.raw.m02; mt.mat.raw.m30 = m.mat.raw.m03;
	mt.mat.raw.m01 = m.mat.raw.m10; mt.mat.raw.m11 = m.mat.raw.m11; mt.mat.raw.m21 = m.mat.raw.m12; mt.mat.raw.m31 = m.mat.raw.m13;
	mt.mat.raw.m02 = m.mat.raw.m20; mt.mat.raw.m12 = m.mat.raw.m21; mt.mat.raw.m22 = m.mat.raw.m22; mt.mat.raw.m32 = m.mat.raw.m23;
	mt.mat.raw.m03 = m.mat.raw.m30; mt.mat.raw.m13 = m.mat.raw.m31; mt.mat.raw.m23 = m.mat.raw.m32; mt.mat.raw.m33 = m.mat.raw.m33;
	return mt;
}

#endif
