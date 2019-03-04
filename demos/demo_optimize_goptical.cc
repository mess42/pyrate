#include <iostream>
#include <fstream>

#include <goptical/core/math/Vector>

#include <goptical/core/material/Base>
#include <goptical/core/material/Conrady>
#include <goptical/core/material/Sellmeier>
#include <goptical/core/material/Abbe>

#include <goptical/core/sys/System>
#include <goptical/core/sys/Lens>
#include <goptical/core/sys/OpticalSurface>
#include <goptical/core/sys/SourcePoint>
#include <goptical/core/sys/Image>

#include <goptical/core/curve/Sphere>
#include <goptical/core/shape/Disk>

#include <goptical/core/trace/Tracer>
#include <goptical/core/trace/Result>
#include <goptical/core/trace/Distribution>
#include <goptical/core/trace/Sequence>
#include <goptical/core/trace/Params>

#include <goptical/core/light/SpectralLine>

#include <goptical/core/analysis/RayFan>
#include <goptical/core/data/Plot>

#include <goptical/core/io/Rgb>
#include <goptical/core/io/RendererSvg>

#include <goptical/core/ref>

using namespace goptical;


int main(void)
{
	sys::system mysys;

	sys::Lens     mylens(math::Vector3(0, 0, 0));

	//       roc,            ap.radius, thickness,
	//
	mylens.add_surface(1/0.031186861,  14.934638, 4.627804137,
        	ref<material::Conrady>::create(1.7, 0.0, 0.0));
	//
	mylens.add_surface(0,              14.934638, 5.417429465);
	//
	mylens.add_surface(1/-0.014065441, 12.766446, 3.728230979,
        	ref<material::Conrady>::create(1.7, 0.0, 0.0));
	//
	mylens.add_surface(1/0.034678487,  11.918098, 4.417903733);
	//
	mylens.add_stop   (                12.066273, 2.288913925);
	mylens.add_surface(0,              12.372318, 1.499288597,
        	ref<material::Conrady>::create(1.7, 0.0, 0.0));
	//
	mylens.add_surface(1/0.035104369,  14.642815, 7.996205852,
        	ref<material::AbbeVd>::create(1.623770, 56.8998));
	//
	mylens.add_surface(1/-0.021187519, 14.642815, 85.243965130);
	//
	mysys.add(mylens);

}
