% Script calculates energy of a single pulse, based on paper:
% D. Giordano, D. Signorino, D. Gallo, H. E. van den Brom, and M. Sira, ‘Methodology for the
% Accurate Measurement of the Power Dissipated by Braking Rheostats’, Sensors, vol. 20, no. 23, Art.
% no. 23, Jan. 2020, doi: 10.3390/s20236935.
% https://www.mdpi.com/1424-8220/20/23/6935
%
% Steps:
% - noise background is fitted
% - fit line is subtracted

% description of indexes of pulses and noise sections:
%     <--> lengthson - durations of on state
%         <--------> lengthsoff - durations of off (noise) state
%               <>        <> tshift_PN
%                 <>    <> tshift_pulse
%      __            __
%     |  |          |  |
%   __|  |__________|  |____
%               ^ ^ ^  ^ ^ ^
%               | | |  | | |
%               | | |  | |  idePN - end of section for noise level estimation
%               | | |  |  ideS - end of pulse shifted by tshift_pulse. also start of section for noise level estimation
%               | | |   ide - end of pulse, really midle of pulse falling slope
%               | |  ids - start of pulse, really middle of pulse rising slope
%               |  idsS - start of pulse shifted by tshift_pulse. also end of section for noise level estimation
%                idsPN - start of section for noise level estimation

function [E1, E2, EPN1, EPN2, Ra1, Ra2, Rb1, Rb2, IaN, IbN] = single_pulse_energy(fs, Vhf, Vdsfhf, Ia, Ib, delta, ideS, idsS, idsPN, idePN, infostr)
        % Fit background %<<<1
        % x axis - just indexes, sufficient for fit:
        x = [idsPN:idsS ideS:idePN];
        % fit noise around pulse for current a
        pa = polyfit(x, [Ia(idsPN:idsS) Ia(ideS:idePN)], 1);
        % tmpIa = Ia(idsS:ideS) - polyval(pa, [idsS:ideS]);
        IaN = polyval(pa, [idsS:ideS]);
        % fit noise around pulse for current b
        pb = polyfit(x, [Ib(idsPN:idsS) Ib(ideS:idePN)], 1);
        % tmpIb = Ib(idsS:ideS) - polyval(pb, [idsS:ideS]);
        IbN = polyval(pb, [idsS:ideS]);

        % resistance of shunts during pulse:
        id = find(max(Ia(idsS:ideS)) == Ia(idsS:ideS));
        id = id(1) + idsS - 1;
        Ra1 = Vhf   (id)./Ia(id);
        Ra2 = Vdsfhf(id)./Ia(id);
        Rb1 = Vhf   (id)./Ib(id);
        Rb2 = Vdsfhf(id)./Ib(id);

        % calculate energy of pulse and noise %<<<1
        % energy of pulse with noise pulse
        E1 = trapz(Vhf   (idsS:ideS).*Ia(idsS:ideS))./fs + trapz(Vdsfhf(idsS:ideS).*Ib(idsS:ideS))./fs;
        E2 = trapz(Vdsfhf(idsS:ideS).*Ia(idsS:ideS))./fs + trapz(Vhf   (idsS:ideS).*Ib(idsS:ideS))./fs;
        % energy of noise part during pulse:
        EPN1 = trapz(Vhf   (idsS:ideS).*IaN)./fs + trapz(Vdsfhf(idsS:ideS).*IbN)./fs;
        EPN2 = trapz(Vdsfhf(idsS:ideS).*IaN)./fs + trapz(Vhf   (idsS:ideS).*IbN)./fs;
        % subtract energy of noise from energy of pulse %<<<1
        E1 = E1 - EPN1;
        E2 = E2 - EPN2;

        % apply pulse shape correction factor K_DC %<<<1
        % load the table:
        %XXX 2DO load from file
        delta_K_DC_table = [0.007 0.02 0.04 0.08 0.010 0.015 0.030; 0.633 0.850 0.923 0.961 0.969 0.979 0.99]';
        % XXX REMOVE! ONLY FOR TESTING:
        delta_K_DC_table = [0 1; 0 1]';
        % interpolate table:
        if delta < min(delta_K_DC_table(:,1)) || delta > max(delta_K_DC_table(:,1))
                error('delta outside K_DC table: delta: %g,  min: %g, max: %g.', ...
                       delta, min(delta_K_DC_table(:, 1)), ...
                       max(delta_K_DC_table(:,1)))
        end
        K_DC = interp1(delta_K_DC_table(:,1), delta_K_DC_table(:, 2), delta);
        % apply factor
        E1 = E1.*K_DC;
        E2 = E2.*K_DC;

        % apply correction for the current sensor K_HOP %<<<1
        % load the table:
        %XXX 2DO load from file
        delta_K_HOP_table = [0.007 0.02 0.08 0.010 0.015 0.030; 1.124 0.960 0.975 0.980 0.987 0.993];
        % XXX REMOVE! ONLY FOR TESTING:
        delta_K_DC_table = [0 1; 0 1]';
        % interpolate table:
        if delta < min(delta_K_HOP_table(:,1)) || delta > max(delta_K_HOP_table(:,1))
                error('delta outside K_HOP table: delta: %g,  min: %g, max: %g.', ...
                       delta, min(delta_K_HOP_table(:, 1)), ...
                       max(delta_K_HOP_table(:,1)))
        end
        K_HOP = interp1(delta_K_HOP_table(:,1), delta_K_HOP_table(:, 2), delta);
        % apply factor
        E1 = E1.*K_HOP;
        E2 = E2.*K_HOP;
endfunction
