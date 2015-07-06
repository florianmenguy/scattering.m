function [phi_ft,energy_sum] = generate_phi_ft(psi_energy_sum,bank_spec)
%%
symmetrized_energy_sum = (psi_energy_sum + psi_energy_sum(end:-1:1)) / 2;
remainder = 1 - min(symmetrized_energy_sum,1);
phi_ft = zeros(size(psi_energy_sum));
signal_dimension = length(bank_spec.size);
original_sizes = bank_spec.size;
if signal_dimension>1
    error('multidimensional phi not ready yet');
end
spin_multiplier = 1 + ~bank_spec.is_spinned;
if bank_spec.is_phi_gaussian
    cutoff_in_dB = 3;
    cutoff = 10^(-cutoff_in_dB/20);
    full_width = bank_spec.phi_bw_multiplier * original_sizes/bank_spec.T;
    frequential_sigma = full_width / sqrt(2*log(1/cutoff));
    range = transpose(0:(bank_spec.size-1));
    numerator = - range.*range;
    denominator = 2 * frequential_sigma^2;
    ln_gaussian = numerator ./ denominator;
    gaussian_decay = exp(ln_gaussian);
    phi_ft = gaussian_decay + gaussian_decay(end:-1:1);
    if bank_spec.is_spinned
        energy_sum = psi_energy_sum + phi_ft;
    else
        energy_sum = psi_energy_sum + 2 * gaussian_decay;
    end
else
    half_support_length = ...
        bank_spec.phi_bw_multiplier/2 * original_sizes/bank_spec.T;
    half_support = 2:half_support_length;
    symmetric_support = original_sizes + 1 - half_support + 1;
    sqrt_truncated_remainder = sqrt(remainder(half_support));
    phi_ft(half_support) = sqrt_truncated_remainder;
    phi_ft(symmetric_support) = sqrt_truncated_remainder;
    energy_sum = psi_energy_sum;
    energy_sum(half_support) = ...
        max(psi_energy_sum(half_support),spin_multiplier);
    if bank_spec.is_spinned
        energy_sum(symmetric_support) = energy_sum(half_support);
    end
end
% We enforce energy conservation for constant signals by setting
% the Fourier transform to exactly 1 at frequency zero.
phi_ft(1+0) = 1;
energy_sum(1+0) = spin_multiplier;
end