module IntegralModule
using LinearAlgebra, SpecialFunctions

# Helper function for overlap
function binomial_factor(factor, o, op, diff_1, diff_2)
    sum = 0.0
    for k in 0:o
        for l in 0:op
            if (factor == (k + l))
                sum = sum + binomial(o,k)*binomial(op,l)*(diff_1)^(o-k) * (diff_2)^(op-l)
            end
        end
    end
    return sum
end

# Overlap integral between primitives
function calOverlapElementGaussian(alpha1, alpha2, o1, o2, pos1, pos2)
    r_average = (alpha1*pos1 + alpha2*pos2)/(alpha1 + alpha2)
    harmonic = (alpha1*alpha2)/(alpha1 + alpha2)
    pos_diff = (pos1[1] - pos2[1])^2 + (pos1[2] - pos2[2])^2 + (pos1[3] - pos2[3])^2
    overlap_fixed = 0
    for m in 0:2:(o1[1] + o2[1])
        for n in 0:2:(o1[2] + o2[2])
            for q in 0:2:(o1[3] + o2[3])
                Cm = binomial_factor(m, o1[1], o2[1], r_average[1] - pos1[1], r_average[1] - pos2[1])
                Cn = binomial_factor(n, o1[2], o2[2], r_average[2] - pos1[2], r_average[2] - pos2[2])
                Cq = binomial_factor(q, o1[3], o2[3], r_average[3] - pos1[3], r_average[3] - pos2[3])
                factorial_m = 1.0
                factorial_n = 1.0
                factorial_q = 1.0
                for s in 0:((m/2)-1)
                    factorial_m = factorial_m*(2*s + 1)
                end
                for t in 0:((n/2)-1)
                    factorial_n = factorial_n*(2*t + 1)
                end
                for u in 0:((q/2)-1)
                    factorial_q = factorial_q*(2*u + 1)
                end
                factorial_binomial_terms = (0.5)^(0.5*(m+n+q)) * factorial_m * factorial_n * factorial_q * Cm * Cn * Cq * pi^1.5 
                alpha_term = (alpha1 + alpha2)^(-1.5 - 0.5*(m+n+q))
                overlap_fixed = overlap_fixed + factorial_binomial_terms * alpha_term * exp(-harmonic*pos_diff)
            end
        end
    end
    return overlap_fixed
end

# Contracted Gaussian Overlap Matrix elements
function calOverlapElementCG(orbital1, orbital2)
    overlap = 0.0
    for i in 1:orbital1.num_gaussians
        for j in 1:orbital2.num_gaussians
            coeff_terms = orbital1.coeffs[i]*orbital2.coeffs[j]*orbital1.norms[i]*orbital2.norms[j]
            overlap_fixed = calOverlapElementGaussian(orbital1.alphas[i], orbital2.alphas[j], orbital1.o, orbital2.o, orbital1.pos, orbital2.pos)
            overlap = overlap + overlap_fixed * coeff_terms 
        end
    end
    return overlap
end

function calKineticEnergyElementGaussian(alpha1, alpha2, o1, o2, pos1, pos2)
    kinetic_term = 0.0
    if (o2[1] >= 2)
        kinetic_term = kinetic_term + o2[1]*(o2[1] - 1)*calOverlapElementGaussian(alpha1, alpha2, o1, [o2[1] - 2, o2[2], o2[3]], pos1, pos2)
    end
    if (o2[2] >= 2)
        kinetic_term = kinetic_term + o2[2]*(o2[2] - 1)*calOverlapElementGaussian(alpha1, alpha2, o1, [o2[1], o2[2] - 2, o2[3]], pos1, pos2)
    end
    if (o2[3] >= 2)
        kinetic_term = kinetic_term + o2[3]*(o2[3] - 1)*calOverlapElementGaussian(alpha1, alpha2, o1, [o2[1], o2[2], o2[3] - 2], pos1, pos2)
    end
    kinetic_term = kinetic_term + 4*(alpha2^2)*(calOverlapElementGaussian(alpha1, alpha2, o1, [o2[1] + 2, o2[2], o2[3]], pos1, pos2) 
    + calOverlapElementGaussian(alpha1, alpha2, o1, [o2[1], o2[2]+2, o2[3]], pos1, pos2) 
    + calOverlapElementGaussian(alpha1, alpha2, o1, [o2[1], o2[2], o2[3]+2], pos1, pos2)) 
    kinetic_term = kinetic_term - 2*alpha2*(2*(o2[1] + o2[2] + o2[3]) + 3)*calOverlapElementGaussian(alpha1, alpha2, o1, o2, pos1, pos2)
    return -0.5*kinetic_term
end

function calKineticEnergyElementCG(orbital1, orbital2)
    kinetic = 0.0
    for i in 1:orbital1.num_gaussians
        for j in 1:orbital2.num_gaussians
            coeff_term = orbital1.coeffs[i]*orbital2.coeffs[j]*orbital1.norms[i]*orbital2.norms[j]
            kinetic_term = calKineticEnergyElementGaussian(orbital1.alphas[i], orbital2.alphas[j], orbital1.o, orbital2.o, orbital1.pos, orbital2.pos)
            kinetic = kinetic + coeff_term*kinetic_term
        end
    end
    return kinetic
end

function CRI(in_x, in_y, in_z, sub_in_x, sub_in_y, sub_in_z, alpha1, alpha2)
    sqrtm1 = sqrt(complex(-1))
    term_1 = (0.5/sqrt(alpha1 + alpha2))^(in_x + in_y + in_z)
    term_2 = sqrtm1^(in_x + in_y + in_z) * factorial(in_x) * factorial(in_y) * factorial(in_z)
    term_3 = (pi/(alpha1 + alpha2))^1.5 * (-1.0)^(sub_in_x + sub_in_y + sub_in_z)
    term_4 = (alpha1 + alpha2)^(-0.5*(in_x + in_y + in_z - 2*(sub_in_x + sub_in_y + sub_in_z)))
    term_denominator = factorial(sub_in_x)*factorial(sub_in_y)*factorial(sub_in_z)
    term_denominator_2 = factorial(in_x - 2*sub_in_x) * factorial(in_y - 2*sub_in_y) * factorial(in_z - 2*sub_in_z)
    cri = (term_1*term_2*term_3*term_4)/(term_denominator*term_denominator_2)
    return cri
end

function Ipk(in_val, sub_in_val, sub_sub_in_val, alpha1, alpha2, pos_diff_component)
    sqrtm1 = sqrt(complex(-1))
    term_1 = (alpha1 + alpha2)^(0.5*(in_val - 2*sub_in_val))
    term_2 = sqrtm1^(in_val - 2*sub_in_val) * factorial(in_val - 2*sub_in_val)
    term_3 = sqrt(pi*(alpha1 + alpha2)) * (-1.0)^(sub_sub_in_val)
    term_4 = (pos_diff_component)^(in_val - 2*sub_in_val - 2*sub_sub_in_val) * (4*(alpha1 + alpha2))^(0.5*(in_val - 2*sub_in_val - 2*sub_sub_in_val))
    term_denominator = factorial(sub_sub_in_val) * factorial(in_val - 2*sub_in_val - 2*sub_sub_in_val)
    ipk = 2*(term_1*term_2*term_3*term_4)/(term_denominator)
    return ipk
end

function Boys(nu, x)
    if (x == 0)
        boys = 1.0/(2*nu + 1)
    else
        p,q = gamma_inc(nu + 0.5, x)
        boys = p*gamma(nu + 0.5)/(2*x^(nu + 0.5))
    end
    return boys
end

function calNuclearAttractionGaussian(alpha1, alpha2, o1, o2, pos1, pos2, posC)
    r_average = (alpha1*pos1 + alpha2*pos2)/(alpha1 + alpha2)
    harmonic = (alpha1*alpha2)/(alpha1 + alpha2)
    pos_diff = (pos1[1] - pos2[1])^2 + (pos1[2] - pos2[2])^2 + (pos1[3] - pos2[3])^2
    posPC = (r_average[1] - posC[1])^2 + (r_average[2] - posC[2])^2 + (r_average[3] - posC[3])^2
    nuclear_attraction_fixed = complex(0.0)
    for k in 0:(o1[1] + o2[1])
        Ck = binomial_factor(k, o1[1], o2[1], r_average[1] - pos1[1], r_average[1] - pos2[1])
        for l in 0:(o1[2] + o2[2])
            Cl = binomial_factor(l, o1[2], o2[2], r_average[2] - pos1[2], r_average[2] - pos2[2])
            for m in 0:(o1[3] + o2[3])
                Cm = binomial_factor(m, o1[3], o2[3], r_average[3] - pos1[3], r_average[3] - pos2[3])
                for q in 0:Int64(floor(k/2))
                    for s in 0:Int64(floor(l/2))
                        for t in 0:Int64(floor(m/2))
                            CRI_val = CRI(k, l, m, q, s, t, alpha1, alpha2)
                            for v in 0:Int64(floor((k - 2*q)/2))
                                Ik_v = Ipk(k, q, v, alpha1, alpha2, r_average[1] - posC[1])
                                for w in 0:Int64(floor((l - 2*s)/2))
                                    Ik_w = Ipk(l, s, w, alpha1, alpha2, r_average[2] - posC[2])
                                    for f in 0:Int64(floor((m - 2*t)/2))
                                        Ik_f = Ipk(m, t, f, alpha1, alpha2, r_average[3] - posC[3])
                                        u_int = Boys(k+l+m - 2*(q + s + t) - v - w - f, (alpha1 + alpha2)*posPC)
                                        nuclear_attraction_fixed = nuclear_attraction_fixed + u_int*Ik_f*Ik_w*Ik_v*CRI_val*Ck*Cl*Cm
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return exp(-harmonic*pos_diff)*nuclear_attraction_fixed/(4*pi^2*(alpha1 + alpha2))
end

function calElectronRepulsionPotential(pos_r, orbital1, orbital2)
    electron_repulsion = complex(0.0)
    for i in 1:orbital1.num_gaussians
        for j in 1:orbital2.num_gaussians
            coeff_term = orbital1.coeffs[i]*orbital2.coeffs[j]*orbital1.norms[i]*orbital2.norms[j]
            electron_value = calNuclearAttractionGaussian(orbital1.alphas[i], orbital2.alphas[j], orbital1.o, orbital2.o, orbital1.pos, orbital2.pos, pos_r)
            electron_repulsion = electron_repulsion + coeff_term*electron_value
        end
    end
    return electron_repulsion
end

function calNuclearAttractionCG(orbital1, orbital2, atoms)
    nuclear_attraction = complex(0.0)
    for c in atoms
        Zc = c.atomic_number
        posC = c.pos 
        for i in 1:orbital1.num_gaussians
            for j in 1:orbital2.num_gaussians
                coeff_term = orbital1.coeffs[i]*orbital2.coeffs[j]*orbital1.norms[i]*orbital2.norms[j]
                nuclear_value = -Zc*calNuclearAttractionGaussian(orbital1.alphas[i], orbital2.alphas[j], orbital1.o, orbital2.o, orbital1.pos, orbital2.pos, posC)
                nuclear_attraction = nuclear_attraction + coeff_term*nuclear_value
            end
        end
    end
    return nuclear_attraction 
end

function theta(factor, o, op, diff_1, diff_2, ri, alphasum)
    cfactor = binomial_factor(factor, o, op, diff_1, diff_2)
    theta = cfactor * factorial(factor) * (alphasum)^(ri - factor)/(factorial(ri)*factorial(factor - 2*ri))
    return theta
end

function gi(lp, lq, rp, rq, f, o1, o2, o3, o4, diff_1, diff_2, diff_3, diff_4, diff_centers, gamma_p, gamma_q)
    theta_terms = (-1.0)^lp * theta(lp, o1, o2, diff_1 ,diff_2, rp, gamma_p) * theta(lq, o3, o4, diff_3, diff_4, rq, gamma_q)
    delta = 0.25*(1.0/gamma_p + 1.0/gamma_q)
    factorial_nume = (-1.0)^f * (2*delta)^(2*(rp + rq)) * factorial(lp + lq - 2*rp - 2*rq) * delta^f * diff_centers^(lp + lq - 2*(rp + rq + f))
    factorial_deno = (4*delta)^(lp + lq) * factorial(f) * factorial(lp + lq - 2*(rp + rq + f))
    gi = theta_terms * factorial_nume/factorial_deno
    return gi
end

function calElectronRepulsionGaussian(alpha1, alpha2, alpha3, alpha4, o1, o2, o3, o4, pos1, pos2, pos3, pos4)
    gamma_p = alpha1 + alpha2
    gamma_q = alpha3 + alpha4
    delta = 0.25*(1.0/gamma_p + 1.0/gamma_q)
    r_average_p = (alpha1*pos1 + alpha2*pos2)/(alpha1 + alpha2)
    r_average_q = (alpha3*pos3 + alpha4*pos4)/(alpha3 + alpha4)
    pos_12 = (pos1[1] - pos2[1])^2 + (pos1[2] - pos2[2])^2 + (pos1[3] - pos2[3])^2
    pos_34 = (pos3[1] - pos4[1])^2 + (pos3[2] - pos4[2])^2 + (pos3[3] - pos4[3])^2
    omega = (2*pi^2/(gamma_p*gamma_q))*sqrt(pi/(gamma_p + gamma_q))*exp(-(alpha1*alpha2*pos_12)/gamma_p)*exp(-(alpha3*alpha4*pos_34)/gamma_q)
    pos_pq = (r_average_p[1] - r_average_q[1])^2 + (r_average_p[2] - r_average_q[2])^2 + (r_average_p[3] - r_average_q[3])^2
    ei_fixed = 0.0
    for lp in 0:(o1[1] + o2[1])
        for rp in 0:Int64(floor(lp/2))
            for lq in 0:(o3[1] + o4[1])
                for rq in 0:Int64(floor(lq/2))
                    for fx in 0:Int64(floor((lp + lq - 2*rp - 2*rq)/2))
                        gx = gi(lp, lq, rp, rq, fx, o1[1], o2[1], o3[1], o4[1], r_average_p[1] - pos1[1], 
                        r_average_p[1] - pos2[1], r_average_q[1] - pos3[1], r_average_q[1] - pos4[1], 
                        r_average_p[1] - r_average_q[1], gamma_p ,gamma_q)
                        for mp in 0:(o1[2] + o2[2])
                            for sp in 0:Int64(floor(mp/2))
                                for mq in 0:(o3[2] + o4[2])
                                    for sq in 0:Int64(floor(mq/2))
                                        for fy in 0:Int64(floor((mp + mq - 2*(sp + sq))/2))
                                            gy = gi(mp, mq, sp, sq, fy, o1[2], o2[2], o3[2], o4[2], r_average_p[2] - pos1[2],
                                            r_average_p[2] - pos2[2], r_average_q[2] - pos3[2], r_average_q[2] - pos4[2],
                                            r_average_p[2] - r_average_q[2], gamma_p, gamma_q)
                                            for np in 0:(o1[3] + o2[3])
                                                for tp in 0:Int64(floor(np/2))
                                                    for nq in 0:(o3[3] + o4[3])
                                                        for tq in 0:Int64(floor(nq/2))
                                                            for fz in 0:Int64(floor((np + nq - 2*(tp + tq))/2))
                                                                gz = gi(np, nq, tp, tq, fz, o1[3], o2[3], o3[3], o4[3], r_average_p[3] - pos1[3],
                                                                r_average_p[3] - pos2[3], r_average_q[3] - pos3[3], r_average_q[3] - pos4[3],
                                                                r_average_p[3] - r_average_q[3], gamma_p, gamma_q)
                                                                boys_term = Boys(lp + lq + mp + mq + np + nq - 2*(rp + rq + sp + sq + tp + tq) - (fx + fy + fz), 
                                                                0.25*pos_pq/delta)
                                                                ei_fixed = ei_fixed + boys_term*gx*gy*gz
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return omega*ei_fixed
end

function calElectronRepulsionCG(orbital1, orbital2, orbital3, orbital4)
    ei = 0.0
    for i in 1:orbital1.num_gaussians
        for j in 1:orbital2.num_gaussians
            for k in 1:orbital3.num_gaussians
                for l in 1:orbital4.num_gaussians
                    coeff_term = orbital1.coeffs[i]*orbital2.coeffs[j]*orbital3.coeffs[k]*orbital4.coeffs[l]
                    norm_term = orbital1.norms[i]*orbital2.norms[j]*orbital3.norms[k]*orbital4.norms[l]
                    G = calElectronRepulsionGaussian(orbital1.alphas[i], orbital2.alphas[j], orbital3.alphas[k], orbital4.alphas[l], 
                    orbital1.o, orbital2.o, orbital3.o, orbital4.o, 
                    orbital1.pos, orbital2.pos, orbital3.pos, orbital4.pos)
                    ei = ei + coeff_term*norm_term*G
                end
            end
        end
    end
    return ei
end

function calIntegrals(total_orbitals, orbitals, Atoms)
    # Initialize Hamiltonian matrices
    overlap_matrix = zeros(Float64, (total_orbitals, total_orbitals))
    kinetic_matrix = zeros(Float64, (total_orbitals, total_orbitals))
    nuclear_matrix = zeros(Float64, (total_orbitals, total_orbitals))
    eri_matrix = zeros(Float64, (total_orbitals, total_orbitals, total_orbitals, total_orbitals))
    
    # Construct Hamiltonian matrices
    for i in 1:total_orbitals
        for j in 1:total_orbitals
            overlap_matrix[i,j] = calOverlapElementCG(orbitals[i], orbitals[j])
            kinetic_matrix[i,j] = calKineticEnergyElementCG(orbitals[i], orbitals[j])
            nuclear_matrix[i,j] = calNuclearAttractionCG(orbitals[i], orbitals[j], Atoms)
            for k in 1:total_orbitals
                for l in 1:total_orbitals
                    eri_matrix[i,k,j,l] = calElectronRepulsionCG(orbitals[i], orbitals[j], orbitals[k], orbitals[l])
                end
            end
        end
    end
    return overlap_matrix, kinetic_matrix, nuclear_matrix, eri_matrix
end

end