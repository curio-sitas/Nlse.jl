function _erk4ip_step!(stepper, model)

	#todo! inplace rk_order + nl fct

	stepper.e = exp.(0.5 * stepper.dz * model.dispersion_term)

	stepper.Uip = stepper.e .* stepper.U

	stepper.k1 = stepper.e .* stepper.NU

	stepper.k2 = model.nonlinear_function(model.ifftp * (stepper.Uip .+ 0.5 * stepper.dz * stepper.k1), model)

	stepper.k3 = model.nonlinear_function(model.ifftp * (stepper.Uip .+ 0.5 * stepper.dz * stepper.k2), model)

	stepper.k4 = model.nonlinear_function(model.ifftp * (stepper.e .* (stepper.Uip .+ stepper.dz * stepper.k3)), model)

	stepper.r = stepper.e .* (stepper.Uip .+ stepper.dz * (stepper.k1 / 6.0 .+ stepper.k2 / 3.0 .+ stepper.k3 / 3.0))

	stepper.U1 = stepper.r .+ stepper.dz * stepper.k4 / 6.0

	stepper.k5 = model.nonlinear_function(model.ifftp * stepper.U1, model)

	stepper.U2 = stepper.r .+ stepper.dz * (stepper.k4 / 15.0 .+ stepper.k5 / 10.0)

end
