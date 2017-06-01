clear
cheboppref.setDefaults('ivpSolver', @chebcolloc2);

% Choose initial R1 and R2
R1 = 0;
R2 = 22;

% Guess the boundary conditions for ur and diff(ur, r)
bc = [0, -1/3];

for k = 1:5
  % Set the interval
  d = [R1, R2];

  % Material parameters
  rho = chebfun(@(r) 2 + cos(2*pi*r/10), d);
  mu  = chebfun(@(r) 5 ./ (1+r/10), d);
  lam  = chebfun(@(r) 10*exp(-(r/10)), d);

  % set up the boundary value problem
  L = chebop(d);
  L.op = @(r, ur) ...
         r.^2 .* rho(r) .* ur  + ...
         r.^2 .* (diff(lam)) .* (diff(ur)) + ...
         2 .* r .* (diff(lam)) .* ur + ...
         r.^2 .* lam(r) .* (diff(ur, 2)) + ...
         2 .* r .* lam(r) .* (diff(ur)) - ...
         2 .* lam(r) .* ur + ...
         2 .* r.^2 .* ((diff(mu))) .* (diff(ur)) + ...
         2 .* r.^2 .* mu(r) .* (diff(ur, 2)) + ...
         4 .* r .* mu(r) .* (diff(ur)) - ...
         4 .* mu(r) .* ur;

  % Specify the boundary condition at R1
  L.lbc = @(ur) [ur-bc(1);diff(ur)-bc(2)];

  % Solve the ODE
  ur = L \ 0;

  % Determine \sigma_{rr} and choose R1 and R2 to be 
  % zeros of this function
  r = chebfun('r', d);
  Srr = lam(r) .* (diff(ur) + ur ./ r + ur ./ r) +...
        2 .* mu(r) .* diff(ur);
  if (k == 1)
    x = roots(Srr);
    if(length(x) < 2)
      error('no roots')
    end
    R1 = x(1);
    R2 = x(2);
    clear x
  end

  % update the boundary conditions
  bc(1) = ur(R1);
  ur_r = diff(ur);
  bc(2) = ur_r(R1);
end

% Display info about version of Chebfun and MATLAB used
fprintf('%% Script run with following version of MATLAB and Chebfun\n');
fprintf(['%% MATLAB Release R' version('-release') '\n'])

chebfid = fopen('Contents.m');
fgetl(chebfid);
str = fgetl(chebfid);
fprintf(['%% Chebfun', str(2:end), '\n']);
fclose(chebfid);

% Check end point stress values
fprintf('\n');
fprintf('%% End point stress = %+e (R1 = %e)\n', Srr(R1), R1);
fprintf('%% End point stress = %+e (R2 = %e)\n', Srr(R2), R2);

% Determine max number of points needed to resolve functions
N = max([length(ur), length(Srr), length(lam), ...
         length(rho), length(mu)]);

% Get the Chebyshev points and barycentric weights
[rn, ~, wn] = chebpts(N, d);

description = {'% Interpolation Points',
               '% Barycentric Weights',
               '% Lame''s first parameter',
               '% Lame''s second parameter / shear modulus',
               '% Density',
               '% radial displacement at t = 0',
               '% radial derivative of radial displacement at t = 0'};
vc = {'rn', 'wn', 'lam', 'mu', 'rho', 'ur', 'ur_r'};
v = [rn, wn, lam(rn), mu(rn), rho(rn), ur(rn), ur_r(rn)];

fprintf('\n');
fprintf('R1 = %23.16e; %% Inner Radius\n', R1);
fprintf('R2 = %23.16e; %% Outer Radius\n', R2);

fprintf('N = %d; %% Number of Chebyshev interpolation points\n', N);

for j = 1:length(vc)
  fprintf('\n%-4s', description{j});
  fprintf('\n%-4s = [', vc{j});
  for k = 1:N
      fprintf('%+23.16e', v(k, j))
    if(mod(k, 2) == 0 && k < N)
      fprintf(';\n        ')
    elseif (k < N)
      fprintf('; ');
    end
  end
  fprintf('];\n')
end

% Set \sigma_{\theta\theta} and \sigma_{\phi\phi} for plotting
Sqq = (2*lam(r) + 2 *mu(r)) .* ur ./ r + lam(r) .* diff(ur);

% Plot the non-zero solution variables
FS = 22;
subplot(1, 2, 1)
plot(chebfun(ur, d), '--', 'LineWidth', 2)
hold on
plot(chebfun(Srr, d), '-', 'LineWidth', 2)
plot(chebfun(Sqq, d), '-.', 'LineWidth', 2)
hold off
xlabel('$r$', 'Interpreter', 'LaTex');
title('Solution', 'Interpreter', 'LaTex')
h_legend = legend('$u_{r}$', '$\sigma_{rr}$', '$\sigma_{qq}$');
set(h_legend, 'FontSize', FS, 'Interpreter', 'LaTex');
set(gca, 'FontSize', FS);

% Plot the material parameters
subplot(1, 2, 2)
plot(chebfun(rho, d), '-', 'LineWidth', 2)
hold on
plot(chebfun(lam, d), '--', 'LineWidth', 2)
plot(chebfun(mu, d), '-.', 'LineWidth', 2)
hold off
h_legend = legend('$\rho$' , '$\lambda$' , '$\mu$');
xlabel('$r$', 'Interpreter', 'LaTex');
title('Material', 'Interpreter', 'LaTex')
set(h_legend, 'FontSize', FS, 'Interpreter', 'LaTex');
set(gca, 'FontSize', FS);
drawnow

