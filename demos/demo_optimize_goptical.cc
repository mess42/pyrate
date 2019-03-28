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
#include <goptical/core/sys/Stop>
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
#include <goptical/core/io/RendererViewport>

#include <goptical/core/ref>

using namespace goptical;


int main(void)
{

    material::Conrady n17(1.7, 0., 0.);
    material::Conrady n15(1.5, 0., 0.);

    shape::Disk lens_shape(2.0); // lens diameter is 2 mm

    sys::system mysys;

    sys::SourcePoint source(sys::SourceAtFiniteDistance, math::Vector3(0., 0., 0.));

    sys::OpticalSurface s1(math::Vector3(0, 0, 2.0),
                           -5.922, 2.0, material::none, n17);
    sys::OpticalSurface s2(math::Vector3(0, 0, 2.0 + 3.0),
                           -3.160, 2.0, n17, material::none);
    sys::OpticalSurface s3(math::Vector3(0, 0, 2.0 + 3.0 + 5.0),
                           15.884, 2.0, material::none, n17);
    sys::OpticalSurface s4(math::Vector3(0, 0, 2.0 + 3.0 + 5.0 + 3.0),
                           -12.756, 2.0, n17, material::none);
    sys::Stop stop(math::Vector3(0, 0, 2.0 + 3.0 + 5.0 + 3.0 + 3.0), 1.0);

    sys::OpticalSurface s6(math::Vector3(0, 0, 2.0 + 3.0 + 5.0 + 3.0 + 3.0 + 2.0),
                           3.125, 2.0, material::none, n15);
    sys::OpticalSurface s7(math::Vector3(0, 0, 2.0 + 3.0 + 5.0 + 3.0 + 3.0 + 2.0 + 3.0),
                           1.479, 2.0, n15, material::none);

    sys::Image image(math::Vector3(0, 0, 2.0 + 3.0 + 5.0 + 3.0 + 3.0 + 2.0 + 3.0 + 19.0), // position
                      10.0);                           // square size,

    mysys.add(source);
    mysys.add(s1);
    mysys.add(s2);
    mysys.add(s3);
    mysys.add(s4);
    mysys.add(stop);
    mysys.add(s6);
    mysys.add(s7);
    mysys.add(image);


    trace::Sequence seq(mysys);
    mysys.get_tracer_params().set_sequential_mode(seq);


    io::RendererSvg renderer("layout.svg", 1024, 100);
    io::RendererViewport &myrenderer = renderer;

    // draw 2d system layout
    mysys.draw_2d_fit(renderer);
    mysys.draw_2d(renderer);

    std::cout << mysys;
    std::cout << seq;

    // trace and draw rays from source

    std::cout << "tracer init" << std::endl;
    trace::tracer tracer(mysys);

    tracer.get_params().set_default_distribution(
        trace::Distribution(trace::MeridionalDist, 21.0));
    tracer.get_trace_result().set_generated_save_state(source);
    tracer.trace();
    tracer.get_trace_result().draw_2d(myrenderer);


}
