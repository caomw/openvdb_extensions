/*
 * Copyright(C) 2016, Blake C. Lucas, Ph.D. (img.science@gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include "SpringlShader.h"
namespace aly {
	void SpringlShader::setForegroundTextureImage(const std::string& textureImage) {
		matcapTextureForeground.load(textureImage, false);
	}
	void SpringlShader::setForegroundTextureImage(const ImageRGBA& textureImage) {
		matcapTextureForeground.load(textureImage, false);
	}
	void SpringlShader::setBackgroundTextureImage(const std::string& textureImage) {
		matcapTextureBackground.load(textureImage, false);
	}
	void SpringlShader::setBackgroundTextureImage(const ImageRGBA& textureImage) {
		matcapTextureBackground.load(textureImage, false);
	}
	SpringlShader::SpringlShader(bool onScreen, const std::shared_ptr<AlloyContext>& context) :
			GLShader(onScreen, context), matcapTextureForeground(onScreen, context) , matcapTextureBackground(onScreen, context) {
		initialize( { },
				R"(
	#version 330
	layout(location = 0) in vec3 vp; 
	layout(location = 1) in vec2 vt; 
	uniform vec4 bounds;
	uniform vec4 viewport;
	out vec2 uv;
	void main() {
	uv=vt;
	vec2 pos=vp.xy*bounds.zw+bounds.xy;
	gl_Position = vec4(2*pos.x/viewport.z-1.0,1.0-2*pos.y/viewport.w,0,1);
	})",
				R"(
	#version 330
	in vec2 uv;
	out vec4 FragColor;
	uniform ivec2 depthBufferSize;
	const float PI=3.1415926535;
	uniform vec4 tint;
	uniform sampler2D matcapTextureForeground;
	uniform sampler2D matcapTextureBackground;
	uniform sampler2D depthForeground;
	uniform sampler2D depthBackground;

	void main() {
	ivec2 pos=ivec2(uv.x*depthBufferSize.x,uv.y*depthBufferSize.y);
	vec4 fg=texelFetch(depthForeground, pos,0);//Do not interpolate depth buffer!
	vec4 bg=texelFetch(depthBackground, pos,0);
	gl_FragDepth=bg.w;
    vec4 rgba;
	if(bg.w<1.0){
      if(abs(fg.w-bg.w)<0.001f){
  	    rgba=tint*texture(matcapTextureForeground,0.5*bg.xy+0.5);    
      } else {
	    rgba=tint*texture(matcapTextureBackground,0.5*bg.xy+0.5);
      }
	} else {
	  rgba=vec4(0.0,0.0,0.0,0.0);
	}
	FragColor=rgba;
	})");

	}
}

