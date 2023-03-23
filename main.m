clc; clear

%%% Wave equation: (d^2/dx^2 +cd d^2/dy^2) p = 1/c^2 (d^2/dt^2) p

%% OR Gate + XOR Gate, Total runtime: 20000
% Domain definition
c = 1;
lambda = 2000; % 10 cm ... 10 cm to 2000 units lambda 
prefactors = (0.1:0.2:1.9); % as a prefactor for frequency sweep
OUTPUT_VALUES = []; % save output pressure values

for i = 1:size(prefactors,2)    
    prefactor = prefactors(i); disp("prefactor is " + num2str(prefactor) + ".");
    freq = prefactor*c/lambda;
    
    tic
    input = 11;
    filename = "OR_input" + num2str(input) + "prefactor_" + num2str(prefactor) + ".gif";
    filename2 = "OR_input" + num2str(input) + "prefactor_" + num2str(prefactor) + ".png";
    out = solver11_OR(lambda, filename, freq, filename2); disp("solver11_OR is completed.");
    OUTPUT_VALUES = [OUTPUT_VALUES; out];
    
    filename = "XOR_input" + num2str(input) + "prefactor_" + num2str(prefactor) + ".gif";
    filename2 = "XOR_input" + num2str(input) + "prefactor_" + num2str(prefactor) + ".png";
    out = solver11_XOR(lambda, filename, freq, filename2); disp("solver11_XOR is completed.");
    OUTPUT_VALUES = [OUTPUT_VALUES; out];
    
   input = 10;
    filename = "OR_input" + num2str(input) + "prefactor_" + num2str(prefactor) + ".gif";
    filename2 = "OR_input" + num2str(input) + "prefactor_" + num2str(prefactor) + ".png";
    out = solver10_OR(lambda, filename, freq, filename2); disp("solver10_OR is completed.");
    OUTPUT_VALUES = [OUTPUT_VALUES; out];
    
    filename = "XOR_input" + num2str(input) + "prefactor_" + num2str(prefactor) + ".gif";
    filename2 = "XOR_input" + num2str(input) + "prefactor_" + num2str(prefactor) + ".png";
    out = solver10_XOR(lambda, filename, freq, filename2); disp("solver10_XOR is completed.");
    OUTPUT_VALUES = [OUTPUT_VALUES; out];
    
     input = 01;
     filename = "OR_input01" + "prefactor_" + num2str(prefactor) + ".gif";
     filename2 = "OR_input01" + "prefactor_" + num2str(prefactor) + ".png";
     out = solver01_OR(lambda, filename, freq, filename2); disp("solver01_OR is completed.");
     OUTPUT_VALUES = [OUTPUT_VALUES; out];
    
    filename = "XOR_input01" + "prefactor_" + num2str(prefactor) + ".gif";
    filename2 = "XOR_input01" + "prefactor_" + num2str(prefactor) + ".png";
    out = solver01_XOR(lambda, filename, freq, filename2); disp("solver01_XOR is completed.");
    OUTPUT_VALUES = [OUTPUT_VALUES; out];
    toc
end

filename = 'savefile.mat';
save(filename)

function out = solver11_OR(lambda, filename, freq, filename2)

    % Geometric constraints  
    l = lambda / 2;
    w = lambda / 10;
    d = 0.1 * lambda;
    d1 = 0.055 * lambda;

    l1 = (0.03 * l);
    l2 = (0.23 * l);
    w1_upper = (0.69 * w); % changing 
    w1_lower = (0.69 * w); % changing
    w3 = (0.05 * w);
    w2_upper = ((w - w1_upper - 4 * w3) / 2);
    w2_lower = ((w - w1_lower - 4 * w3) / 2);
    w4 = ((l-4*l2)/5);
    ld2 = ((l-4*l1)/8);
    ld3 = ((l2-l1)/2);

    Lx = w + l + d1 + d;
    Ly = 2*w + d;
    dx = 1; dy = 1;
    nx = fix(Lx/dx) + 2; ny = fix(Ly/dy) + 2;

    x = linspace(0, Lx, nx); y = linspace(0, Ly, ny);

    % Field definition

    p = 0*ones(nx,ny);
    p_prev = p; p_next = p;

    CFL = 0.7; % c * dt/dx
    c = 1; %c/lambda; freq1 = c/lambda; freq2 = c/lambda;
    dt = CFL*dx/c;
    output_pres = [];
    time = [];
    
    % Initial conditions
    T = 20000*dt; % Runtime 

    t = 0;
    f = figure(1); set(gcf,'color','w');
    set(gcf,'Position',[100 100 1400 500])

    while t<T
        % Absorbing Boundaries
        p_next(1,:) = p(2,:) + ((CFL-1)/(CFL+1))*(p_next(2,:)-p(1,:)); 
        p_next(end,:) = p(end-1,:) + ((CFL-1)/(CFL+1))*(p_next(end-1,:)-p(end,:));

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        t = t+dt;
        p_prev = p; p = p_next; 

        %p_next(1:end,1) = p(1:end,2) + ((CFL-1)/(CFL+1))*(p_next(1:end,2)-p(1:end,1));
        p_next(1:end,end) = p(1:end,end-1) + ((CFL-1)/(CFL+1))*(p_next(1:end,end-1)-p(1:end,end));

        % Source 1
        p(2, end-lambda/10:end-1) = 1*sin(2*pi*t*freq);

        % Source 2
        p(2, 2:lambda/10-1) = 1*sin(2*pi*t*freq);

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        % Finite difference solution
        for i = 2:nx-1 
            for j = 2:ny-1
                p_next(i,j) = 2*p(i,j) - p_prev(i,j) + CFL^2 * ( p(i+1,j) + p(i,j+1) - 4*p(i,j) + p(i-1,j) + p(i,j-1));
            end
        end

        k = abs(p);
        % Rigid Boundaries
        k(1:end, [1 end]) = -1; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        k(1:w+l, w+dx:end-w-2*dx) = -1; % Middle region
        k(w+l+d1:end,[1:w w+d+dx:end]) = -1; % Right bottom and top regions
        k(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = -1; k(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = -1; % 1/3 

        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = -1; % 2/3
        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = -1; % 2/3
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = -1; % 3/3 
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = -1; % 3/3 

        if round(mod(t/dt,25)) == 1
            output = mean(mean(k(end-d/2:end,w+dx:w+d))); %max(max(k(w+l+d1:end,w+dx:w+d)));
            output_pres = [output_pres output];
            time = [time t];
        end
        if round(mod(t/dt,25)) == 1
            % Visualize
            % COLORBAR
            indexValue = 0.4;     % value for which to set a particular color
            topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
            indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
            bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
            largest = max(max(k));
            smallest = min(min(k));
            L = size(p,1);
            index = L*abs(indexValue-smallest)/(largest-smallest);
            customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
                        linspace(bottomcolor(2),indexColor(2),100*index)',...
                        linspace(bottomcolor(3),indexColor(3),100*index)'];
            customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
                        linspace(indexColor(2),topColor(2),100*(L-index))',...
                        linspace(indexColor(3),topColor(3),100*(L-index))'];
            customCMap = [customCMap1;customCMap2];  % Combine colormaps
            colormap(customCMap)

            imagesc(x,y,k'); colorbar; caxis([-1 max(max(abs(k)))]); % caxis([-1 1]); 
            colorbar; axis equal; 
            grid off; box off; axis off    
            ax = gca;
            ax.YDir = 'normal';
            tstring = num2str(t); pstring = num2str(output); astring = num2str(mean(nonzeros(output_pres)));
            txt = "$t = " + tstring + "; P_{out} = " + pstring + "; P_{avg} = " + astring + "$";
            title(txt,'Interpreter','latex'); %title(sprintf('t = %.6f',t));
            drawnow limitrate
            frame = getframe(f); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            if t == dt 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'Loopcount',inf); 
            else 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'WriteMode','append'); 
            end     
        end
    end
    figure(2); set(gcf,'color','w'); plot(time,output_pres); xlabel("Time",'Interpreter','latex'); ylabel("Output Pressure",'Interpreter','latex'); saveas(gcf,filename2);

    out = output_pres;
    close all
end

function out = solver10_OR(lambda, filename, freq, filename2)

    % Geometric constraints  
    l = lambda / 2;
    w = lambda / 10;
    d = 0.1 * lambda;
    d1 = 0.055 * lambda;

    l1 = (0.03 * l);
    l2 = (0.23 * l);
    w1_upper = (0.69 * w); % changing 
    w1_lower = (0.69 * w); % changing
    w3 = (0.05 * w);
    w2_upper = ((w - w1_upper - 4 * w3) / 2);
    w2_lower = ((w - w1_lower - 4 * w3) / 2);
    w4 = ((l-4*l2)/5);
    ld2 = ((l-4*l1)/8);
    ld3 = ((l2-l1)/2);

    Lx = w + l + d1 + d;
    Ly = 2*w + d;
    dx = 1; dy = 1;
    nx = fix(Lx/dx) + 2; ny = fix(Ly/dy) + 2;

    x = linspace(0, Lx, nx); y = linspace(0, Ly, ny);

    % Field definition

    p = 0*ones(nx,ny);
    p_prev = p; p_next = p;

    CFL = 0.7; % c * dt/dx
    c = 1; %c/lambda; freq1 = c/lambda; freq2 = c/lambda;
    dt = CFL*dx/c;
    output_pres = [];
    time = [];
    
    % Initial conditions
    T = 20000*dt; % Runtime 

    t = 0;
    f = figure(1); set(gcf,'color','w');
    set(gcf,'Position',[100 100 1400 500])

    while t<T
        % Absorbing Boundaries
        p_next(1,:) = p(2,:) + ((CFL-1)/(CFL+1))*(p_next(2,:)-p(1,:)); 
        p_next(end,:) = p(end-1,:) + ((CFL-1)/(CFL+1))*(p_next(end-1,:)-p(end,:));

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        t = t+dt;
        p_prev = p; p = p_next; 

        %p_next(1:end,1) = p(1:end,2) + ((CFL-1)/(CFL+1))*(p_next(1:end,2)-p(1:end,1));
        p_next(1:end,end) = p(1:end,end-1) + ((CFL-1)/(CFL+1))*(p_next(1:end,end-1)-p(1:end,end));

        % Source 1
        p(2, end-lambda/10:end-1) = 1*sin(2*pi*t*freq);

        % Source 2
        %p(2, 2:lambda/10-1) = 1*sin(2*pi*t*freq);

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        % Finite difference solution
        for i = 2:nx-1 
            for j = 2:ny-1
                p_next(i,j) = 2*p(i,j) - p_prev(i,j) + CFL^2 * ( p(i+1,j) + p(i,j+1) - 4*p(i,j) + p(i-1,j) + p(i,j-1));
            end
        end

        k = abs(p);
        % Rigid Boundaries
        k(1:end, [1 end]) = -1; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        k(1:w+l, w+dx:end-w-2*dx) = -1; % Middle region
        k(w+l+d1:end,[1:w w+d+dx:end]) = -1; % Right bottom and top regions
        k(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = -1; k(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = -1; % 1/3 

        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = -1; % 2/3
        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = -1; % 2/3
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = -1; % 3/3 
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = -1; % 3/3 

        if round(mod(t/dt,25)) == 1
            output = mean(mean(k(end-d/2:end,w+dx:w+d))); %max(max(k(w+l+d1:end,w+dx:w+d)));
            output_pres = [output_pres output];
            time = [time t];
        end
        if round(mod(t/dt,25)) == 1
            % Visualize
            % COLORBAR
            indexValue = 0.4;     % value for which to set a particular color
            topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
            indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
            bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
            largest = max(max(k));
            smallest = min(min(k));
            L = size(p,1);
            index = L*abs(indexValue-smallest)/(largest-smallest);
            customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
                        linspace(bottomcolor(2),indexColor(2),100*index)',...
                        linspace(bottomcolor(3),indexColor(3),100*index)'];
            customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
                        linspace(indexColor(2),topColor(2),100*(L-index))',...
                        linspace(indexColor(3),topColor(3),100*(L-index))'];
            customCMap = [customCMap1;customCMap2];  % Combine colormaps
            colormap(customCMap)

            imagesc(x,y,k'); colorbar; caxis([-1 max(max(abs(k)))]); % caxis([-1 1]); 
            colorbar; axis equal; 
            grid off; box off; axis off    
            ax = gca;
            ax.YDir = 'normal';
            tstring = num2str(t); pstring = num2str(output); astring = num2str(mean(nonzeros(output_pres)));
            txt = "$t = " + tstring + "; P_{out} = " + pstring + "; P_{avg} = " + astring + "$";
            title(txt,'Interpreter','latex'); %title(sprintf('t = %.6f',t));
            drawnow limitrate
            frame = getframe(f); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            if t == dt 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'Loopcount',inf); 
            else 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'WriteMode','append'); 
            end     
        end
    end
    figure(2); set(gcf,'color','w'); plot(time,output_pres); xlabel("Time",'Interpreter','latex'); ylabel("Output Pressure",'Interpreter','latex'); saveas(gcf,filename2)
    out = output_pres;
    close all
end

function out = solver01_OR(lambda, filename, freq, filename2)

    % Geometric constraints  
    l = lambda / 2;
    w = lambda / 10;
    d = 0.1 * lambda;
    d1 = 0.055 * lambda;

    l1 = (0.03 * l);
    l2 = (0.23 * l);
    w1_upper = (0.69 * w); % changing 
    w1_lower = (0.69 * w); % changing
    w3 = (0.05 * w);
    w2_upper = ((w - w1_upper - 4 * w3) / 2);
    w2_lower = ((w - w1_lower - 4 * w3) / 2);
    w4 = ((l-4*l2)/5);
    ld2 = ((l-4*l1)/8);
    ld3 = ((l2-l1)/2);

    Lx = w + l + d1 + d;
    Ly = 2*w + d;
    dx = 1; dy = 1;
    nx = fix(Lx/dx) + 2; ny = fix(Ly/dy) + 2;

    x = linspace(0, Lx, nx); y = linspace(0, Ly, ny);

    % Field definition

    p = 0*ones(nx,ny);
    p_prev = p; p_next = p;

    CFL = 0.7; % c * dt/dx
    c = 1; %c/lambda; freq1 = c/lambda; freq2 = c/lambda;
    dt = CFL*dx/c;
    output_pres = [];
    time = [];
    
    % Initial conditions
    T = 20000*dt; % Runtime 

    t = 0;
    f = figure(1); set(gcf,'color','w');
    set(gcf,'Position',[100 100 1400 500])

    while t<T
        % Absorbing Boundaries
        p_next(1,:) = p(2,:) + ((CFL-1)/(CFL+1))*(p_next(2,:)-p(1,:)); 
        p_next(end,:) = p(end-1,:) + ((CFL-1)/(CFL+1))*(p_next(end-1,:)-p(end,:));

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        t = t+dt;
        p_prev = p; p = p_next; 

        %p_next(1:end,1) = p(1:end,2) + ((CFL-1)/(CFL+1))*(p_next(1:end,2)-p(1:end,1));
        p_next(1:end,end) = p(1:end,end-1) + ((CFL-1)/(CFL+1))*(p_next(1:end,end-1)-p(1:end,end));

        % Source 1
        %p(2, end-lambda/10:end-1) = 1*sin(2*pi*t*freq);

        % Source 2
        p(2, 2:lambda/10-1) = 1*sin(2*pi*t*freq);

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        % Finite difference solution
        for i = 2:nx-1 
            for j = 2:ny-1
                p_next(i,j) = 2*p(i,j) - p_prev(i,j) + CFL^2 * ( p(i+1,j) + p(i,j+1) - 4*p(i,j) + p(i-1,j) + p(i,j-1));
            end
        end

        k = abs(p);
        % Rigid Boundaries
        k(1:end, [1 end]) = -1; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        k(1:w+l, w+dx:end-w-2*dx) = -1; % Middle region
        k(w+l+d1:end,[1:w w+d+dx:end]) = -1; % Right bottom and top regions
        k(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = -1; k(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = -1; % 1/3 

        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = -1; % 2/3
        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = -1; % 2/3
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = -1; % 3/3 
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = -1; % 3/3 

        if round(mod(t/dt,25)) == 1
            output = mean(mean(k(end-d/2:end,w+dx:w+d))); %max(max(k(w+l+d1:end,w+dx:w+d)));
            output_pres = [output_pres output];
            time = [time t];
        end
        if round(mod(t/dt,25)) == 1
            % Visualize
            % COLORBAR
            indexValue = 0.4;     % value for which to set a particular color
            topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
            indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
            bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
            largest = max(max(k));
            smallest = min(min(k));
            L = size(p,1);
            index = L*abs(indexValue-smallest)/(largest-smallest);
            customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
                        linspace(bottomcolor(2),indexColor(2),100*index)',...
                        linspace(bottomcolor(3),indexColor(3),100*index)'];
            customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
                        linspace(indexColor(2),topColor(2),100*(L-index))',...
                        linspace(indexColor(3),topColor(3),100*(L-index))'];
            customCMap = [customCMap1;customCMap2];  % Combine colormaps
            colormap(customCMap)

            imagesc(x,y,k'); colorbar; caxis([-1 max(max(abs(k)))]); % caxis([-1 1]); 
            colorbar; axis equal; 
            grid off; box off; axis off    
            ax = gca;
            ax.YDir = 'normal';
            tstring = num2str(t); pstring = num2str(output); astring = num2str(mean(nonzeros(output_pres)));
            txt = "$t = " + tstring + "; P_{out} = " + pstring + "; P_{avg} = " + astring + "$";
            title(txt,'Interpreter','latex'); %title(sprintf('t = %.6f',t));
            drawnow limitrate
            frame = getframe(f); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            if t == dt 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'Loopcount',inf); 
            else 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'WriteMode','append'); 
            end     
        end
    end
    figure(2); set(gcf,'color','w'); plot(time,output_pres); xlabel("Time",'Interpreter','latex'); ylabel("Output Pressure",'Interpreter','latex'); saveas(gcf,filename2)
    out = output_pres;
    close all
end

function out = solver11_XOR(lambda, filename, freq, filename2)

    % Geometric constraints  
    l = lambda / 2;
    w = lambda / 10;
    d = 0.1 * lambda;
    d1 = 0.055 * lambda;

    l1 = (0.03 * l);
    l2 = (0.23 * l);
    w1_upper = (0.69 * w); % changing 
    w1_lower = (0.22 * w); % changing
    w3 = (0.05 * w);
    w2_upper = ((w - w1_upper - 4 * w3) / 2);
    w2_lower = ((w - w1_lower - 4 * w3) / 2);
    w4 = ((l-4*l2)/5);
    ld2 = ((l-4*l1)/8);
    ld3 = ((l2-l1)/2);

    Lx = w + l + d1 + d;
    Ly = 2*w + d;
    dx = 1; dy = 1;
    nx = fix(Lx/dx) + 2; ny = fix(Ly/dy) + 2;

    x = linspace(0, Lx, nx); y = linspace(0, Ly, ny);

    % Field definition

    p = 0*ones(nx,ny);
    p_prev = p; p_next = p;

    CFL = 0.7; % c * dt/dx
    c = 1; %c/lambda; freq1 = c/lambda; freq2 = c/lambda;
    dt = CFL*dx/c;
    output_pres = [];
    time = [];
    
    % Initial conditions
    T = 20000*dt; % Runtime 

    t = 0;
    f = figure(1); set(gcf,'color','w');
    set(gcf,'Position',[100 100 1400 500])

    while t<T
        % Absorbing Boundaries
        p_next(1,:) = p(2,:) + ((CFL-1)/(CFL+1))*(p_next(2,:)-p(1,:)); 
        p_next(end,:) = p(end-1,:) + ((CFL-1)/(CFL+1))*(p_next(end-1,:)-p(end,:));

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        t = t+dt;
        p_prev = p; p = p_next; 

        %p_next(1:end,1) = p(1:end,2) + ((CFL-1)/(CFL+1))*(p_next(1:end,2)-p(1:end,1));
        p_next(1:end,end) = p(1:end,end-1) + ((CFL-1)/(CFL+1))*(p_next(1:end,end-1)-p(1:end,end));

        % Source 1
        p(2, end-lambda/10:end-1) = 1*sin(2*pi*t*freq);

        % Source 2
        p(2, 2:lambda/10-1) = 1*sin(2*pi*t*freq);

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        % Finite difference solution
        for i = 2:nx-1 
            for j = 2:ny-1
                p_next(i,j) = 2*p(i,j) - p_prev(i,j) + CFL^2 * ( p(i+1,j) + p(i,j+1) - 4*p(i,j) + p(i-1,j) + p(i,j-1));
            end
        end

        k = abs(p);
        % Rigid Boundaries
        k(1:end, [1 end]) = -1; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        k(1:w+l, w+dx:end-w-2*dx) = -1; % Middle region
        k(w+l+d1:end,[1:w w+d+dx:end]) = -1; % Right bottom and top regions
        k(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = -1; k(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = -1; % 1/3 

        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = -1; % 2/3
        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = -1; % 2/3
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = -1; % 3/3 
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = -1; % 3/3 

        if round(mod(t/dt,25)) == 1
            output = mean(mean(k(end-d/2:end,w+dx:w+d))); %max(max(k(w+l+d1:end,w+dx:w+d)));
            output_pres = [output_pres output];
            time = [time t];
        end
        if round(mod(t/dt,25)) == 1
            % Visualize
            % COLORBAR
            indexValue = 0.4;     % value for which to set a particular color
            topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
            indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
            bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
            largest = max(max(k));
            smallest = min(min(k));
            L = size(p,1);
            index = L*abs(indexValue-smallest)/(largest-smallest);
            customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
                        linspace(bottomcolor(2),indexColor(2),100*index)',...
                        linspace(bottomcolor(3),indexColor(3),100*index)'];
            customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
                        linspace(indexColor(2),topColor(2),100*(L-index))',...
                        linspace(indexColor(3),topColor(3),100*(L-index))'];
            customCMap = [customCMap1;customCMap2];  % Combine colormaps
            colormap(customCMap)

            imagesc(x,y,k'); colorbar; caxis([-1 max(max(abs(k)))]); % caxis([-1 1]); 
            colorbar; axis equal; 
            grid off; box off; axis off    
            ax = gca;
            ax.YDir = 'normal';
            tstring = num2str(t); pstring = num2str(output); astring = num2str(mean(nonzeros(output_pres)));
            txt = "$t = " + tstring + "; P_{out} = " + pstring + "; P_{avg} = " + astring + "$";
            title(txt,'Interpreter','latex'); %title(sprintf('t = %.6f',t));
            drawnow limitrate
            frame = getframe(f); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            if t == dt 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'Loopcount',inf); 
            else 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'WriteMode','append'); 
            end     
        end
    end
    figure(2); set(gcf,'color','w'); plot(time,output_pres); xlabel("Time",'Interpreter','latex'); ylabel("Output Pressure",'Interpreter','latex'); saveas(gcf,filename2)
    out = output_pres;
    close all
end

function out = solver10_XOR(lambda, filename, freq, filename2)

    % Geometric constraints  
    l = lambda / 2;
    w = lambda / 10;
    d = 0.1 * lambda;
    d1 = 0.055 * lambda;

    l1 = (0.03 * l);
    l2 = (0.23 * l);
    w1_upper = (0.69 * w); % changing 
    w1_lower = (0.22 * w); % changing
    w3 = (0.05 * w);
    w2_upper = ((w - w1_upper - 4 * w3) / 2);
    w2_lower = ((w - w1_lower - 4 * w3) / 2);
    w4 = ((l-4*l2)/5);
    ld2 = ((l-4*l1)/8);
    ld3 = ((l2-l1)/2);

    Lx = w + l + d1 + d;
    Ly = 2*w + d;
    dx = 1; dy = 1;
    nx = fix(Lx/dx) + 2; ny = fix(Ly/dy) + 2;

    x = linspace(0, Lx, nx); y = linspace(0, Ly, ny);

    % Field definition

    p = 0*ones(nx,ny);
    p_prev = p; p_next = p;

    CFL = 0.7; % c * dt/dx
    c = 1; %c/lambda; freq1 = c/lambda; freq2 = c/lambda;
    dt = CFL*dx/c;
    output_pres = [];
    time = [];
    
    % Initial conditions
    T = 20000*dt; % Runtime 

    t = 0;
    f = figure(1); set(gcf,'color','w');
    set(gcf,'Position',[100 100 1400 500])

    while t<T
        % Absorbing Boundaries
        p_next(1,:) = p(2,:) + ((CFL-1)/(CFL+1))*(p_next(2,:)-p(1,:)); 
        p_next(end,:) = p(end-1,:) + ((CFL-1)/(CFL+1))*(p_next(end-1,:)-p(end,:));

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        t = t+dt;
        p_prev = p; p = p_next; 

        %p_next(1:end,1) = p(1:end,2) + ((CFL-1)/(CFL+1))*(p_next(1:end,2)-p(1:end,1));
        p_next(1:end,end) = p(1:end,end-1) + ((CFL-1)/(CFL+1))*(p_next(1:end,end-1)-p(1:end,end));

        % Source 1
        p(2, end-lambda/10:end-1) = 1*sin(2*pi*t*freq);

        % Source 2
        %p(2, 2:lambda/10-1) = 1*sin(2*pi*t*freq);

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        % Finite difference solution
        for i = 2:nx-1 
            for j = 2:ny-1
                p_next(i,j) = 2*p(i,j) - p_prev(i,j) + CFL^2 * ( p(i+1,j) + p(i,j+1) - 4*p(i,j) + p(i-1,j) + p(i,j-1));
            end
        end

        k = abs(p);
        % Rigid Boundaries
        k(1:end, [1 end]) = -1; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        k(1:w+l, w+dx:end-w-2*dx) = -1; % Middle region
        k(w+l+d1:end,[1:w w+d+dx:end]) = -1; % Right bottom and top regions
        k(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = -1; k(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = -1; % 1/3 

        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = -1; % 2/3
        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = -1; % 2/3
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = -1; % 3/3 
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = -1; % 3/3 

        if round(mod(t/dt,25)) == 1
            output = mean(mean(k(end-d/2:end,w+dx:w+d))); %max(max(k(w+l+d1:end,w+dx:w+d)));
            output_pres = [output_pres output];
            time = [time t];
        end
        if round(mod(t/dt,25)) == 1
            % Visualize
            % COLORBAR
            indexValue = 0.4;     % value for which to set a particular color
            topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
            indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
            bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
            largest = max(max(k));
            smallest = min(min(k));
            L = size(p,1);
            index = L*abs(indexValue-smallest)/(largest-smallest);
            customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
                        linspace(bottomcolor(2),indexColor(2),100*index)',...
                        linspace(bottomcolor(3),indexColor(3),100*index)'];
            customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
                        linspace(indexColor(2),topColor(2),100*(L-index))',...
                        linspace(indexColor(3),topColor(3),100*(L-index))'];
            customCMap = [customCMap1;customCMap2];  % Combine colormaps
            colormap(customCMap)

            imagesc(x,y,k'); colorbar; caxis([-1 max(max(abs(k)))]); % caxis([-1 1]); 
            colorbar; axis equal; 
            grid off; box off; axis off    
            ax = gca;
            ax.YDir = 'normal';
            tstring = num2str(t); pstring = num2str(output); astring = num2str(mean(nonzeros(output_pres)));
            txt = "$t = " + tstring + "; P_{out} = " + pstring + "; P_{avg} = " + astring + "$";
            title(txt,'Interpreter','latex'); %title(sprintf('t = %.6f',t));
            drawnow limitrate
            frame = getframe(f); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            if t == dt 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'Loopcount',inf); 
            else 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'WriteMode','append'); 
            end     
        end
    end
    figure(2); set(gcf,'color','w'); plot(time,output_pres); xlabel("Time",'Interpreter','latex'); ylabel("Output Pressure",'Interpreter','latex'); saveas(gcf,filename2)
    out = output_pres;
    close all
end

function out = solver01_XOR(lambda, filename, freq, filename2)

    % Geometric constraints  
    l = lambda / 2;
    w = lambda / 10;
    d = 0.1 * lambda;
    d1 = 0.055 * lambda;

    l1 = (0.03 * l);
    l2 = (0.23 * l);
    w1_upper = (0.69 * w); % changing 
    w1_lower = (0.22 * w); % changing
    w3 = (0.05 * w);
    w2_upper = ((w - w1_upper - 4 * w3) / 2);
    w2_lower = ((w - w1_lower - 4 * w3) / 2);
    w4 = ((l-4*l2)/5);
    ld2 = ((l-4*l1)/8);
    ld3 = ((l2-l1)/2);

    Lx = w + l + d1 + d;
    Ly = 2*w + d;
    dx = 1; dy = 1;
    nx = fix(Lx/dx) + 2; ny = fix(Ly/dy) + 2;

    x = linspace(0, Lx, nx); y = linspace(0, Ly, ny);

    % Field definition

    p = 0*ones(nx,ny);
    p_prev = p; p_next = p;

    CFL = 0.7; % c * dt/dx
    c = 1; %c/lambda; freq1 = c/lambda; freq2 = c/lambda;
    dt = CFL*dx/c;
    output_pres = [];
    time = [];
    
    % Initial conditions
    T = 20000*dt; % Runtime 

    t = 0;
    f = figure(1); set(gcf,'color','w');
    set(gcf,'Position',[100 100 1400 500])

    while t<T
        % Absorbing Boundaries
        p_next(1,:) = p(2,:) + ((CFL-1)/(CFL+1))*(p_next(2,:)-p(1,:)); 
        p_next(end,:) = p(end-1,:) + ((CFL-1)/(CFL+1))*(p_next(end-1,:)-p(end,:));

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        t = t+dt;
        p_prev = p; p = p_next; 

        %p_next(1:end,1) = p(1:end,2) + ((CFL-1)/(CFL+1))*(p_next(1:end,2)-p(1:end,1));
        p_next(1:end,end) = p(1:end,end-1) + ((CFL-1)/(CFL+1))*(p_next(1:end,end-1)-p(1:end,end));

        % Source 1
        %p(2, end-lambda/10:end-1) = 1*sin(2*pi*t*freq);

        % Source 2
        p(2, 2:lambda/10-1) = 1*sin(2*pi*t*freq);

        % Rigid Boundaries
        p(1:end, [1 end]) = 0; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        p(1:w+l, w+dx:end-w-2*dx) = 0; % Middle region
        p(w+l+d1:end,[1:w w+d+dx:end]) = 0; % Right bottom and top regions
        p(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = 0; p(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = 0; % 1/3 

        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = 0; % 2/3
        p([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = 0; % 2/3
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = 0; % 3/3 
        p([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = 0; % 3/3 

        % Finite difference solution
        for i = 2:nx-1 
            for j = 2:ny-1
                p_next(i,j) = 2*p(i,j) - p_prev(i,j) + CFL^2 * ( p(i+1,j) + p(i,j+1) - 4*p(i,j) + p(i-1,j) + p(i,j-1));
            end
        end

        k = abs(p);
        % Rigid Boundaries
        k(1:end, [1 end]) = -1; % Upper and lower boundaries
        %p([1 end], 1:end) = 0; % Left and right boundaries
        k(1:w+l, w+dx:end-w-2*dx) = -1; % Middle region
        k(w+l+d1:end,[1:w w+d+dx:end]) = -1; % Right bottom and top regions
        k(w:w+l, [1+dx:w3+dx w-w3-dx:w]) = -1; k(w:w+l, [end-w3:end-dx end-w-dx:end-w+w3-dx]) = -1; % 1/3 

        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [w3+2*dx:w3+dx+w2_upper w-w3-1*dx-w2_upper:w-w3-2*dx]) = -1; % 2/3
        k([w:w+w4 w+w4+l2:w+2*w4+l2 w+2*w4+2*l2:w+3*w4+2*l2 w+3*w4+3*l2:w+4*w4+3*l2 w+4*w4+4*l2:w+5*w4+4*l2], [end-w3-w2_lower:end-w3-dx end-w+w3+0*dx:end-w+w3+w2_lower+1*dx]) = -1; % 2/3
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [2*dx+w2_upper+w3:w3+2*dx+w2_upper+w3 w-w3-w2_upper-w3-2*dx:w-w2_upper-w3-2*dx]) = -1; % 3/3 
        k([w:w+w4+ld3 w+w4+l2-ld3:w+2*w4+l2+ld3 w+2*w4+2*l2-ld3:w+3*w4+2*l2+ld3 w+3*w4+3*l2-ld3:w+4*w4+3*l2+ld3 w+4*w4+4*l2-ld3:w+5*w4+4*l2], [end-w3-w2_lower-w3-1*dx:end-w2_lower-w3-1*dx end-w+2*dx+w2_lower+w3:end-w+w3+2*dx+w2_lower+w3]) = -1; % 3/3 

        if round(mod(t/dt,25)) == 1
            output = mean(mean(k(end-d/2:end,w+dx:w+d))); %max(max(k(w+l+d1:end,w+dx:w+d)));
            output_pres = [output_pres output];
            time = [time t];
        end
        if round(mod(t/dt,25)) == 1
            % Visualize
            % COLORBAR
            indexValue = 0.4;     % value for which to set a particular color
            topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
            indexColor = [1 1 1];       % color for indexed data value (white = [1 1 1])
            bottomcolor = [0 0 1];      % color for minimum data value (blue = [0 0 1])
            largest = max(max(k));
            smallest = min(min(k));
            L = size(p,1);
            index = L*abs(indexValue-smallest)/(largest-smallest);
            customCMap1 = [linspace(bottomcolor(1),indexColor(1),100*index)',...
                        linspace(bottomcolor(2),indexColor(2),100*index)',...
                        linspace(bottomcolor(3),indexColor(3),100*index)'];
            customCMap2 = [linspace(indexColor(1),topColor(1),100*(L-index))',...
                        linspace(indexColor(2),topColor(2),100*(L-index))',...
                        linspace(indexColor(3),topColor(3),100*(L-index))'];
            customCMap = [customCMap1;customCMap2];  % Combine colormaps
            colormap(customCMap)

            imagesc(x,y,k'); colorbar; caxis([-1 max(max(abs(k)))]); % caxis([-1 1]); 
            colorbar; axis equal; 
            grid off; box off; axis off    
            ax = gca;
            ax.YDir = 'normal';
            tstring = num2str(t); pstring = num2str(output); astring = num2str(mean(nonzeros(output_pres)));
            txt = "$t = " + tstring + "; P_{out} = " + pstring + "; P_{avg} = " + astring + "$";
            title(txt,'Interpreter','latex'); %title(sprintf('t = %.6f',t));
            drawnow limitrate
            frame = getframe(f); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            if t == dt 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'Loopcount',inf); 
            else 
                imwrite(imind,cm,filename,'gif','DelayTime',0.01, 'WriteMode','append'); 
            end     
        end
    end
    figure(2); set(gcf,'color','w'); plot(time,output_pres); xlabel("Time"); ylabel("Output Pressure"); saveas(gcf,filename2)
    out = output_pres;
    close all
end