import { useEffect, useMemo, useRef, useState, useLayoutEffect } from "react";
import * as d3 from "d3";

const MARGIN = { top: 30, right: 30, bottom: 50, left: 50 };


const useDimensions = (targetRef) => {
    const getDimensions = () => {
        return {
            width: targetRef.current ? targetRef.current.offsetWidth : 0,
            height: targetRef.current ? targetRef.current.offsetHeight : 0
        };
    };

    const [dimensions, setDimensions] = useState(getDimensions);

    const handleResize = () => {
        setDimensions(getDimensions());
    };

    useEffect(() => {
        window.addEventListener("resize", handleResize);
        return () => window.removeEventListener("resize", handleResize);
    }, []);

    useLayoutEffect(() => {
        handleResize();
    }, []);

    return dimensions;
}


const getBarplotData = (featureValues, selectedPoints, selectedGroupName) => {
    const uniqueValues = [...new Set(featureValues)];
    var selectedCounts = {};
    var otherCounts = {};
    uniqueValues.forEach((value) => {
        selectedCounts[value] = 0;
        otherCounts[value] = 0;
    });

    var isSelected = new Array(featureValues.length).fill(false);
    selectedPoints.forEach((point) => {
        isSelected[point] = true;
    });

    for (var i = 0; i < featureValues.length; i++) {
        var value = featureValues[i];
        if (isSelected[i]) {
            selectedCounts[value] = selectedCounts[value] + 1;
        } else {
            otherCounts[value] = otherCounts[value] + 1;
        }
    }

    // normalize by dividing by the number of selected/other points (selectedPoints.length)
    uniqueValues.forEach((value) => {
        selectedCounts[value] = selectedCounts[value] / (selectedPoints.length == 0 ? 1 : selectedPoints.length);
        otherCounts[value] = otherCounts[value] / (featureValues.length - selectedPoints.length);
    });

    selectedCounts["x"] = selectedGroupName;
    otherCounts["x"] = "other";

    return [selectedCounts, otherCounts];
}

export const StackedBarplot = ({
    featureValues,
    xlabel,
    selectedPoints,
    selectedGroupName,
    colorMap
}) => {
    // bounds = area inside the graph axis = calculated by substracting the margins
    const axesRef = useRef(null);
    const targetRef = useRef(null);

    const dimensions = useDimensions(targetRef);
    const width = dimensions.width;
    const height = dimensions.height;
    const boundsWidth = width - MARGIN.right - MARGIN.left;
    const boundsHeight = height - MARGIN.top - MARGIN.bottom;


    const data = getBarplotData(featureValues, selectedPoints, selectedGroupName);
    const allGroups = data.map((d) => String(d.x));

    const stackSeries = d3.stack().keys(colorMap["ticks"]);
    const series = stackSeries(data);
    console.log(series);
    console.log(`colormap: ${colorMap["colors"]}`);


    // Y axis
    const yScale = useMemo(() => {
        return d3
            .scaleLinear()
            .domain([0, 1])
            .range([boundsHeight, 0]);
    }, [data, height]);

    // X axis
    const xScale = useMemo(() => {
        return d3
            .scaleBand()
            .domain(allGroups)
            .range([0, boundsWidth])
            .padding(0.05);
    }, [data, width]);

    // Color Scale
    var colorScale = d3
        .scaleOrdinal()
        .domain(colorMap["ticks"])
        .range(colorMap["colors"]);

    // Render the X and Y axis using d3.js, not react
    useEffect(() => {
        const svgElement = d3.select(axesRef.current);
        svgElement.selectAll("*").remove();
        const xAxisGenerator = d3.axisBottom(xScale);
        svgElement
            .append("g")
            .attr("transform", "translate(0," + boundsHeight + ")")
            .call(xAxisGenerator);

        const yAxisGenerator = d3.axisLeft(yScale);
        svgElement.append("g").call(yAxisGenerator);

        svgElement.append("text")
            .attr("y", boundsHeight + 40)
            .attr("x", boundsWidth / 2)
            .style("text-anchor", "middle")
            .text(xlabel);
    }, [xScale, yScale, boundsHeight]);

    const rectangles = series.map((subgroup, i) => {
        return (
            <g key={i}>
                {subgroup.map((group, j) => {
                    return (
                        <rect
                            key={j}
                            x={xScale(group.data.x)}
                            y={yScale(group[1])}
                            height={(yScale(group[0]) - yScale(group[1])) < 0 ? 0 : yScale(group[0]) - yScale(group[1])}
                            width={xScale.bandwidth()}
                            fill={colorScale(subgroup.key)}
                            opacity={0.9}
                        ></rect>
                    );
                })}
            </g>
        );
    });

    return (
        <div className='flex flex-wrap items-center justify-left text-sm min-w-full max-w-full h-[240px]'
            ref={targetRef}>
            <svg width={width} height={height}>
                <g
                    width={boundsWidth}
                    height={boundsHeight}
                    transform={`translate(${[MARGIN.left, MARGIN.top].join(",")})`}
                >
                    {rectangles}
                </g>
                <g
                    width={boundsWidth}
                    height={boundsHeight}
                    ref={axesRef}
                    transform={`translate(${[MARGIN.left, MARGIN.top].join(",")})`}
                />
            </svg>
        </div>
    );
};
