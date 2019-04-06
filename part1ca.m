%Part 1b: Takes into account correlated noise

%Top function.
function[D_PRIME NEURONS] = part1ca()
    
    for j=1:300
        f_max(j) = gamrnd(20/0.6, 0.6);
    end
    
    [D_PRIME offset] = return_d_prime_vs_numberofneurons(pi/3, pi/2, 1, f_max, 100, 2, 300)
    
    %for K
    %for i=1:951
    %    K(i) = 0.5 + (i-1)*0.01;
    %end
    
    %for NEURONS
    for i=1:300-offset
        NEURONS(i) = i+offset;
    end
    
 
  
end

%Function that returns pop vector to given stimulus parameters phi and c
%(representing orientation and contrast respectively)
function [pop_vec] = von_mis_neurons(phi, c, neurons, f_max, K)
    f_0=5;    %offset aka baseline
    
    phi_step=pi/(neurons-1); %increments by which preferred orientations between neurons differ
    pop_vec = 1:neurons; %to store the population vector
    for i=0:neurons-1
        phi_i = i*phi_step; %preferred orientation of each neuron
        fm = f_max(i+1);
        f_each = f_0 + c*fm*exp(K*(cos(2*(phi-phi_i)) - 1)); %output of each neuron
        %f_each_with_noise = poissrnd(f_each); %output with noise added - not done yet
        pop_vec(i+1) = f_each; %add the output to the population vector
    end
   
end

function[variance_matrix] = return_variance_matrix(neurons)
    variance_matrix = zeros(neurons, neurons);
    for i=1:neurons
        for j=1:neurons
            variance_matrix(i, j) = 5; %because f_0 is 5 for all neurons            
        end
    end
end

function[correlation_matrix] = return_correlation_matrix(neurons)
    correlation_matrix = zeros(neurons, neurons);
    
    c_max = 0.3;
    tau = 0.5;
    
    for i=1:neurons
        for j=1:neurons
            phi_i = (i-1)*pi/(neurons-1);
            phi_j = (j-1)*pi/(neurons-1);
            correlation_matrix(i, j) = c_max*exp(-abs(phi_i - phi_j)/tau);
        end
    end
end

function [cov_matrix] = return_cov_matrix(phi1, phi2, neurons)
    cov_matrix = zeros(neurons, neurons);
    
    correlation_matrix = return_correlation_matrix(neurons);
    variance_matrix = return_variance_matrix(neurons);
    
    cov_matrix = correlation_matrix.*variance_matrix;
    
end

function [w] = weight_calculation(phi1, phi2, c, neurons, f_max, K)
    f1 = von_mis_neurons(phi1, c, neurons, f_max, K);
    f2 = von_mis_neurons(phi2, c, neurons, f_max, K);
    
    %Covariance matrix calculation
    C = return_cov_matrix(phi1, phi2, neurons); 
    C_inverse = inv(C);
    
    %Add correlated Gaussian noise
    f1 = mvnrnd(f1, C);
    f2 = mvnrnd(f2, C);
    
    delta_f = f1-f2;


    w = C_inverse*(delta_f)';
end

%Function that calculates weight vectors and preferred orientations for multiple
%values of delta_phi. Returns the weight vectors, phiprefs and delta_phis
function[weight_vectors phipref delta_phi] = weight_vs_phiPref(phi1, c, neurons, number_of_trials, K)
    
    for i=1:neurons
        phipref(i)=(i-1)*pi/(neurons-1);
    end
    
    for m=1:number_of_trials
        phi2 = (m-1)*pi/(number_of_trials);
        w = weight_calculation(phi1, phi2, c, neurons, K);
        for q=1:neurons
            weight_vectors(m,q) = w(q);
        end
        delta_phi(m) = phi2-phi1;
    end
    

end


%Function that calculates weight vectors and preferred orientations for two
%values of delta_phi. Returns the weight vectors, phiprefs and delta_phis
function[weight_vectors phipref delta_phi] = weight_vs_phiPref_reduced()
    phi1 = 0;
    phi2 = pi/2;
    c=1;
    
    for i=1:100
        phipref(i)=(i-1)*pi/99;
    end
    
    w = weight_calculation(phi1, phi2, c);
    
    for q=1:100
        weight_vectors(1,q) = w(q);
    end
    
    delta_phi(1) = phi2-phi1;    

end

%Function that plots weight vectors vs. preferred orientations for multiple
%values of delta_phi.
function[] = plot_weight_vs_phiPref(weight_vectors, phipref)
    for m=1:size(weight_vectors,1)
        plot(phipref, weight_vectors(m,1:100))
        hold on
    end
end

%Function that plots weight vector values as color plot for preferred
%orientation and delta_phi. Unlike the above function which takes as input
%the output of the weight_vs_phiPref function (a list of weight_vectors and
%a list of preferred orientations),
%this function takes phi1, c, number_of_neurons, number_of_trials as input,
%calls the weight_vs_phiPref function multiple times, and plots the results
function[] = plot_weight_vs_deltaPhi_vs_phiPref(phi1, c, neurons, number_of_trials, K)
    [weight_vectors phi_prefs delta_phis] = weight_vs_phiPref(phi1, c, neurons, number_of_trials, K);
    pcolor(phi_prefs, delta_phis, weight_vectors);
end

%Plots the magnitude of w vs. deltaPhi
function[] = weight_vs_deltaPhi()
    [W phiPref dPHI] = weight_vs_phiPref();
    
    for m=1:size(W,1)
        w_mags(m) = norm(W(m,1:100));
    end
    
    plot(dPHI, w_mags);
    
end

%Calculates the decision variable for a pair of input orientations using
%one of the orientations for 
function [d] = return_decision_variable(phi1, phi2, c, chosen_phi, neurons, f_max, K)
    w = weight_calculation(phi1, phi2, c, neurons, f_max, K);
    
    f1 = von_mis_neurons(chosen_phi, c, neurons, f_max, K);
    
    C = return_cov_matrix(phi1, phi2, neurons);
    f1 = mvnrnd(f1, C);
    
    d = dot(w, f1);

end

%Gives decision variables for many trials
function [D] = decision_variable_multiple(phi1, phi2, c, neurons, f_max, K, number_of_trials)

    for i=1:number_of_trials
        D(1,i) = return_decision_variable(phi1, phi2, c, phi1, neurons, f_max, K);
        D(2,i) = return_decision_variable(phi1, phi2, c, phi2, neurons, f_max, K);
    end
end

%Calculates d_prime from decision variables across many trials
function [d_prime] = return_d_prime(phi1,phi2,c, neurons, f_max, K, number_of_trials)
    D = decision_variable_multiple(phi1, phi2, c, neurons, f_max, K, number_of_trials);
    
    d1mean = mean(D(1,1:size(D,2)));
    d2mean = mean(D(2,1:size(D,2)));
    
    d1var = var(D(1,1:size(D,2)));
    d2var = var(D(2,1:size(D,2)));
    
    d_prime = (d1mean - d2mean)/(sqrt(d1var + d2var));
    
end

%Function that returns multiple d_primes for different numbers of neurons,
%keeping K constant
function [D_PRIME offset] = return_d_prime_vs_numberofneurons(phi1, phi2, c, f_max, number_of_trials, neurons_lower, neurons_upper)

    offset=neurons_lower-1;
    
    for i=neurons_lower:neurons_upper        
        D_PRIME(i-offset) = return_d_prime(phi1, phi2, c, i, f_max, 1, number_of_trials); %Keep K constant at 1
        disp('run no.');
        disp(i-offset);
    end
end

%Function that returns multiple d_primes for different K, keeping number
%of neurons constant
function [D_PRIME] = return_d_prime_vs_K(phi1, phi2, c, number_of_trials, K_lower, K_upper)

    s=0.01;
    for i=1:951
        D_PRIME(i) = return_d_prime(phi1, phi2, c, 100, 0.5 + (i-1)*s, number_of_trials); %Keep number of neurons constant at 100
        disp('run no.');
        disp(i);
    end
end