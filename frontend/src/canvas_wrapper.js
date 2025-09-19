import { useState, useEffect } from "react";
import { Canvas } from "./canvas"

export default function CanvasWrapper({ setScatterLoaded, setScatterplot }) {
    const [isLoading, setIsLoading] = useState(true)

    useEffect(() => {
        setIsLoading(false)
    }, [])

    if (isLoading) {
        return <p>...</p>;
    }

    return (
        <Canvas
            setScatterplot={setScatterplot}
            setScatterLoaded={setScatterLoaded}
        />
    )
}