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

#ifndef INCLUDE_SPRINGLSHADER_H_
#define INCLUDE_SPRINGLSHADER_H_

#include <CommonShaders.h>
namespace aly {
	class SpringlShader: public GLShader {
	private:
		GLTextureRGBA matcapTextureForeground;
		GLTextureRGBA matcapTextureBackground;

	public:
		SpringlShader( bool onScreen = true, const std::shared_ptr<AlloyContext>& context = AlloyDefaultContext());
		void setForegroundTextureImage(const std::string& textureImage);
		void setForegroundTextureImage(const ImageRGBA& textureImage);
		void setBackgroundTextureImage(const std::string& textureImage);
		void setBackgroundTextureImage(const ImageRGBA& textureImage);

		template<class T, int C, ImageType I> void draw(const GLTexture<T, C, I>& imageForeground,const GLTexture<T, C, I>& imageBackground, CameraParameters& camera, const box2px& bounds,
				const box2px& viewport, const RGBAf& tint = RGBAf(1, 1, 1, 1)) {
			begin() .set("matcapTextureForeground", matcapTextureForeground, 0)
					.set("matcapTextureBackground", matcapTextureBackground, 1)
					.set("depthForeground", imageForeground, 2)
					.set("depthBackground", imageBackground, 3)
					.set("depthBufferSize",  imageForeground.dimensions())
					.set("bounds",bounds)
					.set("viewport", viewport)
					.set("tint", tint).draw( imageForeground).end();
		}
	};

}
#endif /* INCLUDE_SPRINGLSHADER_H_ */
