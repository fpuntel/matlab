
# Inicio
clear -all
clock

# Variaveis principais
n = 2
p = 5
randomSupim = 0.5
randomObs = 0.01
raio = 15;
totalExperimentos = 50;

# Matriz Verdade
x=linspace(0,2*pi,n);
y=linspace(0,2*pi,n);
for ind = 1:n
  for ind2 = 1:n
    xt(ind, ind2) = sin(x(ind))*sin(y(ind2)); 
  endfor
endfor

# Construcao Matriz Verdade (xt)
indLinear = 1;
for ind = 1:n
  for ind2 = 1:n
    xt_linear(indLinear) = xt(ind, ind2);
    indLinear = indLinear + 1;
  endfor
endfor

# Construcao Matriz Background (xb)
indLinear = 1;
for ind = 1:n
  for ind2 = 1:n
    xb(ind, ind2) = normrnd (xt(ind,ind2), randomSupim);
    xb_linear(indLinear) = xb(ind, ind2);
    indLinear = indLinear + 1;
  endfor
endfor
xb_linear = xb_linear';

# Localização dos Vetores de observacao, observacao_x e observacao_y
for indObs = 1:p
  obsX(indObs) = rand(1) * n;
  if(obsX(indObs) < 1)
    obsX(indObs) = obsX(indObs) + 1;
  endif
  obsY(indObs) = rand(1) * n;
  if(obsY(indObs) < 1)
    obsY(indObs) = obsY(indObs) + 1;
  endif
endfor

#clear obsValue;
#for indObs = 1:p
#  oX = int32(obsX(indObs));
#  oY = int32(obsY(indObs));
  
#  obsValue(indObs) = normrnd (xt(oX,oY), randomObs);

  # Somente para conferencia
#  obsTrue(indObs) = xt(oX, oY);
#  obsBackGround(indObs) = xb(oX, oY);
 
  # obsLinear posicao linear da observacao 
  # utilizado no calculo de R
#  obsLinear(indObs) = ((oX-1) * n) + oY;
  #xt(oX,oY)
  #xt_linear(obsLinear(indObs))
#endfor

# Matriz H
indColuna = 1; # Pois a matriz H deve ser pX(n*n)

for indObs = 1:p
  Sum(indObs) = 0;
endfor

for indObs = 1:p
#disp('xxxxxx')
  for indX = 1:n
  #disp('x')
    for indY = 1:n
      distancia = ((obsX(indObs) - indX) ^ 2) + ((obsY(indObs) - indY) ^ 2);
      distancia = sqrt(distancia);
      #distancia
      if distancia == 0
        H(indObs, indColuna) = 1e300;      
      else
        H(indObs, indColuna) = (1/distancia)^10;
      endif  
      Sum(indObs) = Sum(indObs) + H(indObs, indColuna);
      indColuna = indColuna + 1;
      
    endfor 
  endfor
  indColuna = 1;
endfor

for indObs = 1:p
  soma = 0;
  for ind = 1:n*n
    H(indObs, ind) = H(indObs, ind) / Sum(indObs);
  endfor
endfor

disp("Tamanho matriz H");
size(H)

# Matriz B
for indErros = 1:totalExperimentos
  
  for ind = 1:n*n
    temp_xb_linear(ind) =  normrnd (xt_linear(ind), randomSupim);
  endfor
  
  for ind = 1:n*n
    diferencaEB(ind, indErros) = temp_xb_linear(ind) - xt_linear(ind);
  endfor
endfor

for ind = 1:n*n # Matriz media
  aux = 0;
  for indErros = 1:totalExperimentos
    aux = aux + diferencaEB(ind, indErros); # Neste caso varia o indice das linhas, pois é feito a media de todos os itens
  endfor
  mediaErros(ind) = aux/totalExperimentos;
endfor
mediaErros = mediaErros';

clear B;
B = 0;
for indErros = 1:totalExperimentos # Construcao Matriz B
  vetor = diferencaEB(:, indErros);
  aux = (vetor - mediaErros) * ((vetor - mediaErros)');
  B = B + aux;
endfor
B = B / totalExperimentos;

disp("\Tamanho matriz B")
size(B)

# Construcao Matriz R

# Definiçao dos valores de observaçao
truth_nas_Obs = H * xt_linear';
disp("\H")
size(H)
H
disp("\xt_linear")
size(xt_linear')
xt_linear'
for indObs = 1:p
  obsValue(indObs) = normrnd (truth_nas_Obs(indObs), randomObs);
endfor
  
  
for indErros = 1:totalExperimentos
  
  for indObs = 1:p
    temp_obs_linear(indObs) = normrnd (truth_nas_Obs(indObs), randomObs);
  endfor
 
  for indObs = 1:p
    diferencaER(indObs, indErros) = temp_obs_linear(indObs) - truth_nas_Obs(indObs);
  endfor  
endfor

diferencaER

clear mediaErros;
mediaErros = 0;
for indObs = 1:p # Media
  aux = 0;
  for indErros = 1:totalExperimentos
    aux = aux + diferencaER(indObs, indErros); # Soma todas as colunas
  endfor
  mediaErros(indObs) = aux/totalExperimentos; 
endfor
mediaErros = mediaErros';

clear R;
R = 0;
for indErros = 1:totalExperimentos
  vetor = diferencaER(:, indErros);
  aux = (vetor - mediaErros) * ((vetor - mediaErros)');
  R = R + aux;
endfor
R = R/totalExperimentos;

disp("\Tamanho matriz R")
size(R)

R

# Calculos Finais 

HBH_T = H * B * (H');
BH_T = B * (H');
Hxb = H * xb_linear;

# Blue
K = HBH_T + R;
K = BH_T * inv(K);
inovacao = K * (obsValue' - Hxb);
size(inovacao)
xa_linear = xb_linear + inovacao;

indLinear = 1;
for ind = 1:n
  for ind2 = 1:n
    xa(ind, ind2) = xa_linear(indLinear);
    indLinear = indLinear + 1;
  endfor
endfor

# Visualizações...
figure(1)
hold off
surf(xt)
hold on
plot3(obsX,obsY,obsValue,'*')

figure(2)
surf(xb)

figure(3)
surf(xa)


#### Para conferencia dos resultados finais
ind = 1;
xt_linear = xt_linear';

erroGlobal = 0;
for ind = 1:n
  for ind2 = 1:n
	  erroGlobal = erroGlobal + (xt(ind,ind2) - xb(ind,ind2) ) ^2;
  endfor
endfor
disp("\nErro background")
erroGlobal/(n*n)

erroGlobal = 0;
for ind = 1:n
  for ind2 = 1:n
	  erroGlobal = erroGlobal + (xt(ind,ind2) - xa(ind,ind2) ) ^2;
  endfor
endfor
disp("\nErro analise")
erroGlobal/(n*n)

erroGlobal = 0;
for a =1: p
	erroGlobal = erroGlobal + (truth_nas_Obs(a) - obsValue(a) ) ^2;
endfor
disp("\nerro obs")
erroGlobal/p

erroGlobal = 0;
#xa_linear = xa_linear';
xa_nas_Obs = H * xa_linear;
for a =1: p	
	erroGlobal = erroGlobal + (truth_nas_Obs(a) - xa_nas_Obs(a) ) ^2;
endfor
disp("\nerro analise nos pontos de observação")
erroGlobal/p

clock
