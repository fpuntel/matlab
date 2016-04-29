# Limpeza dados e horario
clear -all
clock

# Variaveis principais
n = 1500;                   # Tamanho do background
p = 400;                   # Total de observacoes
randomSupim = 0.5;        # Random para background
randomObs = 0.03;          # Random para as observacoes
totalExperimentos = 500;  # Total experimentos para ErrorB e ErrorR

# Vetores
xb = 0;
xt_linear = 0;
K = 0;

# Matriz Verdade
x=linspace(0,2*pi,n);
xt_linear=sin(x); 
size(xt_linear);

# Background
for ind = 1:n
  xb(ind) = normrnd (xt_linear(ind), randomSupim);
endfor

# Observacoes
for indObs = 1:p
  
  posObs(indObs) = (rand(1) * n); # n*n pois eh um vetor
  if posObs(indObs) < 1
    posObs(indObs) = posObs(indObs) + 1;
  endif
  
  posRealObs = uint8(posObs(indObs));
  obsValue(indObs) = normrnd (xt_linear(posRealObs), randomObs);
  
  # Somente para conferencia
  obsTrue(indObs) = xt_linear(posRealObs);
  obsBackGround(indObs) = xb(posRealObs);
  
endfor

# Matriz H
for c =1: p
	Sum(c) = 0;
endfor

for indObs = 1:p # OBS
  for ind = 1:n  # Background
    distancia = (posObs(indObs) - ind)^2;
    distancia = sqrt(distancia);
   if distancia == 0
      H(indObs, ind) = 1e300;
      #  Colocado um valor alto, para caso, a obs coincidir com a posicao 
			#  exata do ponto na matriz, assim após dividir cada item pela soma 
			#  de todos, o ponto sofrerá uma maior influência dos que os demais.
    else
      H(indObs, ind) = (1/distancia)^12;
    endif
    Sum(indObs) = Sum(indObs) + H(indObs, ind);
  endfor  
endfor

for indObs = 1:p
  for ind = 1: n
   H(indObs, ind) = H(indObs, ind) / Sum(indObs);
  endfor
endfor

# Matriz B
for indErros = 1:totalExperimentos

  for ind = 1: n
    tempo_xb(ind) = normrnd (xt_linear(ind), randomSupim);
  endfor
  
  for ind = 1:n
    diferencaEB(ind, indErros) = tempo_xb(ind) - xt_linear(ind) ;
    # Neste caso varia o indice de colunas. Pois cada linha representa um teste
  endfor
  
endfor

for ind = 1:n # Media dos experimentos
  aux = 0;
  for indErros = 1:totalExperimentos
    aux = aux + diferencaEB(ind, indErros);
    # Neste caso varia o indice das linhas, pois é feito a media de todos os itens
  endfor
  mediaErros(ind) = aux/totalExperimentos;
endfor
mediaErros = mediaErros';

clear B; 
B=0;
for indErros = 1:totalExperimentos # Construcao matriz B
  vetor = diferencaEB(:, indErros);
  aux = (vetor - mediaErros) * ((vetor - mediaErros)');
  B = B + aux;
endfor
B = B / totalExperimentos;

disp("Matriz B")
size(B)

# Construcao Matriz R
truth_nas_Obs = H * xt_linear';
for indErros = 1:totalExperimentos
  
  for indObs = 1:p
    temp_obs_linear(indObs) = normrnd (xt_linear(uint8(posObs(indObs))), randomObs);
  endfor
  
  for indObs = 1:p
    diferencaER(indObs, indErros) = temp_obs_linear(indObs) - truth_nas_Obs(indObs);
  endfor
endfor

for indObs = 1:p # Construcao Matriz Media
  aux = 0;
  for indErros = 1:totalExperimentos
    aux = aux + diferencaER(indObs, indErros); # Soma todas as colunas
  endfor
  mediaErros(indObs) = aux/totalExperimentos;
endfor
mediaErros = mediaErros';

clear R;
R=0;
for indErros = 1:totalExperimentos
  vetor = diferencaER(:, indErros);
  aux = (vetor - mediaErros) * ((vetor - mediaErros)');
  R = R + aux;
endfor
R = R/totalExperimentos;

disp("Matriz R")
size(R)

# Calculos Finais
HBH_T = H * B * (H'); #Matriz pxp
BH_T = B * (H');
xb = xb';
Hxb = H * xb;

K = HBH_T + R;
K = BH_T * inv(K);
inovacao = K * (obsValue' - Hxb);
xa = xb + inovacao;

hold off
figure(1)
plot(xt_linear)
hold on
plot(xb, 'g')
hold on
plot(xa, 'm')
hold on
for indObs = 1:p
  
  plot(posObs(indObs),obsValue(indObs),'*');
  
endfor


#### Para conferencia dos resultados finais
erroGlobal = 0;
for ind = 1:n
	erroGlobal = erroGlobal + (xt_linear(ind) - xb(ind) ) ^2;
endfor
disp("\nErro background")
erroGlobal/n

erroGlobal = 0;
for ind = 1:n
	erroGlobal = erroGlobal + (xt_linear(ind) - xa(ind) ) ^2;
endfor
disp("\nErro analise")
erroGlobal/n



