global EPSILON   = 1e-3
global NUM_ITER  = 1e5

function newton_basins (f, df, l, u, p)
  # Acessa as variaveis globais
  global EPSILON;
  global NUM_ITER;

  # Remove os valores da dupla
  [ l1, l2 ] = unpack_double(l);
  [ u1, u2 ] = unpack_double(u);
  [ p1, p2 ] = unpack_double(p);

  # Calcula a densidade do pixel em relação aos eixos real e imaginário
  deltaHorizontal = (1.0 * (l2 - l1) / (p1 - 1));
  deltaVertical = (1.0 * (u2 - u1) / (p2 - 1));

  # Cria um vetor que armazena todas as raízes encontradas
  roots = zeros(1e5, 1);
  n = 0;

  # Defini o arquivo de saida
  filename = "output.txt";
  output = fopen (filename, "w");

  for h = 0:(p1 - 1);
    for v = 0:(p2 - 1);
      # Mapeia o dominio em pixels
      x0 = (l1 + h * deltaHorizontal) + i * (u1 + v * deltaVertical);

      # Calcula o método de newton nesse pixel
      [ success, x ] = newton (f, df, x0, EPSILON, NUM_ITER);

      if success
        rootIndex = -1;

        # Procura, se existir, a raiz correspondente de x
        for j = 1:n
          root = roots(j);

          if (abs(x - root) <= 2 * EPSILON)
            rootIndex = j;
            break
          endif
        endfor

        # Se não encontrou, armazena a sua existência
        if rootIndex == -1
          n = n + 1;
          roots(n) = x;
          rootIndex = n;
        endif

        # Escreve no arquivo de saída a raiz
        write (h, v, rootIndex, output);

      else
        # Escreve no arquivo de saída que não convergiu (raiz 0)
        write (h, v, 0, output);
      endif

      # Para testes apenas
      % printf("Evaluating for (%f + %fi): %d - (%f + %fi)\n",
      %     real(x0), imag(x0), success, real(x), imag(x));
    endfor

    # Para testes apenas
    # break
  endfor

  # Fecha o arquivo
  fclose (output);

endfunction

function [ successful, x ] = newton (f, df, x0, epsilon, num_iter)
  x = x0;
  successful = false;

  while (num_iter > 0)
    # Atualiza o contador de iterações
    num_iter = num_iter - 1;

    # Calcula f(x) e f'(x)
    fx = f(x);
    dfx = df(x);

    # Se |f(x)| <= eps, então x0 convergiu para x
    if (abs(fx) <= epsilon)
      successful = true;
      break
    endif

    # Se |f'(x)| <= eps, então x0 não converge para nenhuma raiz
    if (abs(dfx) <= epsilon)
      break
    endif

    # Note que não precisamos verificar que dfx é diferente de 0,
    # porque se for, ele cai no caso de cima (não converge)
    x = x - fx/dfx;
  endwhile
endfunction

# Imprime a raiz
function write (x, y, root, output)
  fprintf (output, "%d %d %d\n", x, y, root);
endfunction


# Útil para desempacotar duplas (como Octave não tem isso na linguagem?)
function [ x0, x1 ] = unpack_double (x)
  x0 = x(1);
  x1 = x(2);
endfunction


# Uso:
l = [ -2 2 ];
u = [ -2 2 ];
p = [ 1024 768 ];
f = inline ("x^4 - 1");
df = inline ("4*(x^3)");
newton_basins (f, df, l, u, p)
