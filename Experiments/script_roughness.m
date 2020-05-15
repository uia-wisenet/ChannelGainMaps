% Roughness experiment
n_x = 10;
x_span = 10;
v_x = sort(x_span*rand(n_x, 1));
v_y = sort(rand(n_x, 1));

m_Dr = diag([v_x(1); diff(v_x)]);
m_Ls  = zeros(n_x); % for the second derivative
for ii = 1:n_x
    m_Ls(ii:end, ii) = m_Dr(ii,ii);
end

m_Lf = m_Dr/2; % for the first derivative
for ii = 1:n_x-1
    m_Lf(ii+1:end, ii) = m_Lf(ii+1, ii+1)+m_Lf(ii, ii);
end

m_Lsi = inv(m_Ls);
m_Lfi = inv(m_Lf);

m_M = (m_Lfi'*m_Lsi')*m_Dr*(m_Lsi*m_Lfi);
roughness = 1/2*v_y'*m_M*v_y;


n_samples = 10000;
vus_x = linspace(0, x_span, n_samples); %uniform samples of x
stepsize = mean(diff(vus_x));
v_yf = m_Lf\v_y; %first derivative
v_ys = m_Ls\v_yf;  %second derivative
vus_ys = zeros(size(vus_x));
for ii = n_x:-1:1
    vus_ys(vus_x<=v_x(ii))=v_ys(ii);
end
vus_yf = cumsum(stepsize*vus_ys);
vus_y  = cumsum(stepsize*vus_yf);

figure(1001); clf
stem(v_x, v_y);
hold on
plot(vus_x, vus_y)

figure(1002); clf
plot(vus_x, vus_yf); hold on
plot(v_x, v_yf);

figure(1003); clf
plot(vus_x, vus_ys); hold on
stem(v_x, v_ys);