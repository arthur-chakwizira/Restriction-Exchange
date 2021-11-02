function resex_mc_spin_visualiser(film)
%This function visualises particle positions at every step
R_new = film.R_new;
states_new = film.states_new;
save_film = film.save_film;

set(film.plot_1, 'XData', R_new(states_new==1,1)*1e6, 'YData', R_new(states_new==1,2)*1e6);
set(film.plot_2, 'XData', R_new(states_new==2,1)*1e6, 'YData', R_new(states_new==2,2)*1e6);
drawnow

if save_film
    writerObj = film.writerObj;
    frame = getframe;
    writeVideo(writerObj,frame);
end
end