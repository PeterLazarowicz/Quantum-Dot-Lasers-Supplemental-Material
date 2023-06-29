function z_dot = TPM_1F_d_eq(t,z,p)
% Returns the time derivative of the Two Particle Model with no
% two-fermion correlations and different quantum dots (TPM_1F_d) using complex
% variables. The variable z passed to this function contains either the
% real and imaginary parts of the complex fields, or are the amplitude of
% the real fields.  The mapping is as follows:
%
% CIM variables [Mark's thesis equations (3.76)-(3.80) and equations (1-5)
% of Francesco's SlowlyVaryingEqs.pdf note.
%
% b = z(1)+1i*z(2);
% Bb = z(3);
% Cc = z(5+(1:n));
% Vc = z(5+n+(1:n))+1i*z(5+2*n+(1:n));
% bCv = z(5+3*n+(1:n))+1i*z(5+4*n+(1:n)); 
%
% Boson-fermion and boson-boson operators [Mark's thesis equations (3.81),
% (3.82) and (3.83) and equations (6-8) of Francesco's SlowlyVaryingEqs.pdf
% note. 
% bCc = z(5+5*n+(1:n))+1i*z(5+6*n+(1:n));
% bVc = z(5+7*n+(1:n))+1i*z(5+8*n+(1:n));
% bb = z(4)+1i*z(5);
%
% where we have used the notation b = <b> and B = <b†>, etc.
%
% The structure p contains the equation coefficients

persistent t_write old_pump
if (isempty(t_write))
  t_write = 0.0;
  old_pump = p.r;
end
if (p.r ~= old_pump)
  t_write = 0;
  old_pump = p.r;
end

% Number of quantum dots 
n = p.N;

% CIM variables [Mark's thesis equations (3.76)-(3.80) and equations (1-5)
% of Francesco's SlowlyVaryingEqs.pdf note.
b = z(1)+1i*z(2);
Bb = z(3);
Cc = z(5+(1:n));
Vc = z(5+n+(1:n))+1i*z(5+2*n+(1:n));
bCv = z(5+3*n+(1:n))+1i*z(5+4*n+(1:n)); 

% Boson-fermion and boson-boson operators [Mark's thesis equations (3.81),
% (3.82) and (3.83) and equations (6-8) of Francesco's SlowlyVaryingEqs.pdf
% note. 
bCc = z(5+5*n+(1:n))+1i*z(5+6*n+(1:n));
bVc = z(5+7*n+(1:n))+1i*z(5+8*n+(1:n));
bb = z(4)+1i*z(5);

% The TPM equations for the expectation values [Francesco's notes
% 2F_EV_7-8.pdf]
%  + 2*real(sum(p.g.*bCv))

b_dot = -p.gamma_c*b + sum(conj(p.g).*Vc);

Bb_dot = -2*p.gamma_c*Bb + 2*real(sum(p.g.*bCv));

Vc_dot = -(p.gamma-1i*p.Delta_nu).*Vc + p.g.*(2*bCc - b);

Cc_dot = p.r*(1-Cc) - (p.gamma_nr + p.gamma_nl)*Cc - 2*real(p.g.*bCv);

bCv_dot = -(p.gamma + p.gamma_c + 1i*p.Delta_nu).*bCv ...
  + conj(p.g).*(Cc + Bb.*(2*Cc-1) + 4*real(b.*conj(bCc)) - 4.*b.*conj(b).*Cc) ...
  + conj((Vc*(Vc')-diag(Vc*(Vc')))*p.g);

bCc_dot =  -(p.gamma_nr + p.gamma_c)*bCc - b*Cc*p.gamma_nl  ...
  - p.g.*(conj(Vc).*bb + 2*b.*bCv - 2*b^2.*conj(Vc)) ...
  - conj(p.g).*(Vc.*Bb + b.*conj(bCv) - 2*abs(b)^2.*Vc ...
  + conj(b).*bVc) + p.r*b.*(1-Cc) ...
  + (Cc*Vc.'-diag(Cc*Vc.'))*conj(p.g);

bVc_dot = -(p.gamma + p.gamma_c-1i*p.Delta_nu).*bVc ...
  + p.g.*(bb.*(2*Cc-1) + 4*b.*bCc - 4*Cc.*b^2) ...
  + (Vc*Vc.'-diag(Vc*Vc.'))*conj(p.g);

bb_dot = -2*p.gamma_c*bb + 2*p.g'*bVc;

% Package all the time derivatives in a single vector
z_dot = [real(b_dot); imag(b_dot); ...
  Bb_dot; ...
  real(bb_dot); imag(bb_dot); ...
  Cc_dot; ...
  real(Vc_dot); imag(Vc_dot); ...
  real(bCv_dot); imag(bCv_dot); ...
  real(bCc_dot); imag(bCc_dot); ...
  real(bVc_dot); imag(bVc_dot)];
if t <=10
    load('ev.mat','Bb_ev')
    print =b_dot;
    Bb_ev = [Bb_ev print];
    save('ev.mat',"Bb_ev")
end
