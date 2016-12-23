#include <iostream>
#include <fstream>

#include <Goptical/Math/Vector>

#include <Goptical/Material/Base>
#include <Goptical/Material/Sellmeier>

#include <Goptical/Sys/System>
#include <Goptical/Sys/OpticalSurface>
#include <Goptical/Sys/SourcePoint>
#include <Goptical/Sys/Image>

#include <Goptical/Curve/Sphere>
#include <Goptical/Shape/Disk>

#include <Goptical/Trace/Tracer>
#include <Goptical/Trace/Result>
#include <Goptical/Trace/Distribution>
#include <Goptical/Trace/Sequence>
#include <Goptical/Trace/Params>

#include <Goptical/Light/SpectralLine>

#include <Goptical/Analysis/RayFan>
#include <Goptical/Data/Plot>

#include <Goptical/Io/Rgb>
#include <Goptical/Io/RendererSvg>

using namespace Goptical;


int main(void)
{
	Sys::System sys

	// code from examples/tessar_lens/tessar.cc:70
	Sys::Lens     lens(Math::Vector3(0, 0, 0));

	//       roc,            ap.radius, thickness,
	//
	lens.add_surface(1/0.031186861,  14.934638, 4.627804137,
        	ref<Material::Conrady>::create(1.7, 0.0, 0.0));
	//
	lens.add_surface(0,              14.934638, 5.417429465);
	//
	lens.add_surface(1/-0.014065441, 12.766446, 3.728230979,
        	ref<Material::Conrady>::create(1.7, 0.0, 0.0));
	//
	lens.add_surface(1/0.034678487,  11.918098, 4.417903733);
	//
	lens.add_stop   (                12.066273, 2.288913925);
	lens.add_surface(0,              12.372318, 1.499288597,
        	ref<Material::Conrady>::create(1.7, 0.0, 0.0));
	//
	lens.add_surface(1/0.035104369,  14.642815, 7.996205852,
        	ref<Material::AbbeVd>::create(1.623770, 56.8998));
	//
	lens.add_surface(1/-0.021187519, 14.642815, 85.243965130);
	//
	sys.add(lens);

}
