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

	sys::SourcePoint source(sys::SourceAtInfinity, math::Vector3(0., 0., 1.));

	sys::Lens     mylens(math::Vector3(0, 0, 2.0));

	//       roc,            ap.radius, thickness,
	//
	mylens.add_surface(-5.922,  2.0, 3.0,
        	ref<material::Conrady>::create(1.7, 0.0, 0.0));
	//
	mylens.add_surface(-3.160,  2.0, 5.0, material::none);
	//
	mylens.add_surface(15.884, 2.0, 3.0,
        	ref<material::Conrady>::create(1.7, 0.0, 0.0));
	//
	mylens.add_surface(-12.756,  2.0, 3.0, material::none);
	//
	mylens.add_stop   (                1.0, 2.0);
	mylens.add_surface(3.125,              2.0, 3.0,
        	ref<material::Conrady>::create(1.5, 0.0, 0.0));
	//
	mylens.add_surface(1.479,  2.0, 19.0, material::none);
	//
	mylens.add_surface(0.0, 10.0, 0.0, material::none);
	//
	mysys.add(mylens);

	trace::Sequence seq(mysys);
	mysys.get_tracer_params().set_sequential_mode(seq);


	io::RendererSvg renderer("layout.svg", 1024, 100);

	// draw 2d system layout
	mysys.draw_2d_fit(renderer);
	mysys.draw_2d(renderer);

	std::cout << mysys;
	std::cout << seq;

	// trace and draw rays from source

	std::cout << "tracer init" << std::endl;
	trace::tracer tracer(mysys);
	std::cout << "1" << std::endl;

	tracer.get_params().set_default_distribution(
		trace::Distribution(trace::MeridionalDist, 1.0));
	tracer.get_trace_result().set_generated_save_state(source);
	tracer.trace();
	tracer.get_trace_result().draw_2d(renderer);


}
