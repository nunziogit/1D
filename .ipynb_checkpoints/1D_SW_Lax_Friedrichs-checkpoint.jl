using Plots, Plots.Measures
default(size = (1200, 400), framestyle = :box, label = false, grid = false, margin = 10mm, lw = 6, labelfontsize = 20, tickfontsize = 20, titlefontsize = 24)

# Include the external file with the plotting function
include("plot_results.jl")


@views function Lax_Friedrichs_1D()
	# physics
	lx      = 30.0
	gravit  = 9.81
	hinitl  = 1.0
	uinitl  = 0.0
	hinitr  = 0.5
	uinitr  = 0.0
	timeout = 1.2
	# numerics
	nx   = 100
	nt   = 10000
	#nvis = 20
	cfl  = 0.9
	# derived numerics
	dx = lx / nx
	xc = LinRange(dx / 2, lx - dx / 2, nx)
	igate = nx ÷ 2
	# array initialisation
	h = zeros(Float64, nx)
	u = zeros(Float64, nx)
	flux = zeros(Float64, 2, nx - 1)
	flux_diff = zeros(Float64, nx - 2)  # Preallocate for diff results
	fl = zeros(Float64, 2)
	fr = zeros(Float64, 2)
	for i ∈ 1:nx
		#h[i] , u[i]) =  i < igate ? (hinitl, uinitl) : (hinitr, uinitr)
		h[i] = i < igate ? hinitl : hinitr
		u[i] = i < igate ? uinitl : uinitr
	end
	#@views h,u

	q = @. u * h
	h_i = copy(h)
	u_i = copy(u)
	# time loop
	time = 0.0
	timeout_reached = false
	for it ∈ 1:nt
		# set CFL condition
		λ1 = abs.(@. u .- sqrt(gravit * h))
		λ2 = abs.(@. u .+ sqrt(gravit * h))
		λmax = maximum(max.(λ1, λ2))

		#λ[1, :] = @. u - sqrt(gravit * h)
		#λ[2, :] = @. u + sqrt(gravit * h)
		#λmax = maximum(abs.(λ))
		dt = cfl * dx / λmax
		if time + dt > timeout
			dt = timeout - time
			timeout_reached = true
		end
		time += dt
		for i ∈ 1:nx-1
			hl = h[i]
			ul = u[i]
			ql = hl * ul
			hr = h[i+1]
			ur = u[i+1]
			qr = hr * ur
			fl[1] = hl * ul
			fl[2] = ul^2 * hl + 0.5 * gravit * hl^2
			fr[1] = hr * ur
			fr[2] = ur^2 * hr + 0.5 * gravit * hr^2
			flux[1, i] = 0.5 * (fl[1] + fr[1]) - 0.5 * dx / dt * (hr - hl)
			flux[2, i] = 0.5 * (fl[2] + fr[2]) - 0.5 * dx / dt * (qr - ql)
		end
		#UPDATE

		h[2:end-1] .-= dt ./ dx .* (diff(flux[1, :]))
		q[2:end-1] .-= dt ./ dx .* (diff(flux[2, :]))
		u = @. q / h

		# Exit the loop if timeout_reached is true
		if timeout_reached
			println("Timeout reached")
			break
		end
	end
	# Call the plot_results function after the loop
    plot_results(xc, h_i, h, u_i, u, lx, timeout)
end

Lax_Friedrichs_1D()
