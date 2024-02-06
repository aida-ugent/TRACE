'use client'

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
        <Canvas data={{
            "filename": "none",
            "data": {
                "x": [0, 0.02, 0.01],
                "y": [0, 0, 0.02],
                "z": [0, 1, 2]
            }
        }} setScatterplot={setScatterplot} setScatterLoaded={setScatterLoaded} />
    )
}