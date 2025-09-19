import { useCallback, useRef, useEffect } from 'react';
import createScatterplot from 'regl-scatterplot';
export var scatterplot;

// need canvas to draw the scatterplot, onMount advice from https://stackoverflow.com/questions/33924150/how-to-access-canvas-context-in-react
export function Canvas({ setScatterLoaded, setScatterplot}) {
    const canvasRef = useRef(null);

    const refCallback = useCallback((canvas) => {
        if (canvas == null) return;
        canvasRef.current = canvas;

        if (!scatterplot) {
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
                actionKeyMap: { remove: 'ctrl', merge: 'cmd', lasso: 'shift' },
                showReticle: true,
                reticleColor: [.3, .3, .3, 0.66],
                cameraDistance: 1.2,
                backgroundColor: '#ffffff',
                pointColor: ["#EBAC23", "#B80058", "#008CF9"],
                antiAliasing: 0.8,
            });
        }
    }, []);

    // Handle data drawing separately
    useEffect(() => {
        // Canvas initialization only - data is handled in scatterplot.js
        if (scatterplot && canvasRef.current) {
            setScatterLoaded(true);
            setScatterplot(scatterplot);
        }
    }, []);

    return (
        <div className="relative w-full h-full">
            <canvas
                id="canvas"
                ref={refCallback}
                className={`bg-white w-full h-full`}
            />
        </div>
    );
};