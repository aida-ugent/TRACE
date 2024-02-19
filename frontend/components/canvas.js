'use-client'

import { useCallback } from 'react';
import createScatterplot from 'regl-scatterplot';
export var scatterplot;

// need canvas to draw the scatterplot, onMount advice from https://stackoverflow.com/questions/33924150/how-to-access-canvas-context-in-react
export function Canvas({ data, setScatterLoaded, setScatterplot}) {
    const refCallback = useCallback((canvas) => {
        if (canvas == null) return;
        //const { width, height } = canvas.getBoundingClientRect();
        scatterplot = createScatterplot({
            canvas,
            width: 'auto',
            height: 'auto',
            pointSize: 5,
        });

        scatterplot.set({
            colorBy: 'z',
            opacityBy: 'density',
            pointSizeSelected: 3,
            opacityInactiveScale: 0.9,
            pointOutlineWidth: 3,
            lassoColor: [.3, .3, .3, .7],
            lassoMinDelay: 0,
            lassoMinDist: 1,
            keyMap: { ctrl: 'merge', shift: 'lasso', alt: 'rotate' },
            showReticle: true,
            reticleColor: [.3, .3, .3, 0.66],
            cameraDistance: 1.2,
            backgroundColor: '#ffffff',
            pointColor: ["#EBAC23", "#B80058", "#008CF9"]
        });

        scatterplot.draw({
            x: data.data['x'],
            y: data.data['y'],
            z: data.data['z']
        })
        .then(() => {
            setScatterLoaded(true);
            setScatterplot(scatterplot);
        });
    }, [])

    return (

        <canvas
            ref={refCallback}
            className='bg-white'
        />

    )
}


