import React, { useState, useEffect, useMemo, useRef } from 'react';

export const RangeFilter = ({ 
  featureValues, 
  title, 
  onRangeChange, 
  minValue,
  maxValue
}) => {
  const [isDragging, setIsDragging] = useState(null); // 'min', 'max', or null
  const sliderRef = useRef(null);

  // Calculate the actual min and max from the feature values
  const { actualMin, actualMax } = useMemo(() => {
    if (!featureValues || featureValues.length === 0) {
      return { actualMin: 0, actualMax: 1 };
    }
    const min = Math.min(...featureValues);
    const max = Math.max(...featureValues);
    return { actualMin: min, actualMax: max };
  }, [featureValues]);

  const getPercentage = (value) => {
    if (actualMax === actualMin) return 0;
    return ((value - actualMin) / (actualMax - actualMin)) * 100;
  };

  const getValueFromPercentage = (percentage) => {
    return actualMin + (percentage / 100) * (actualMax - actualMin);
  };

  const handleMouseDown = (handle) => (e) => {
    e.preventDefault();
    setIsDragging(handle);
  };

  const handleMouseMove = (e) => {
    if (!isDragging || !sliderRef.current) return;

    const rect = sliderRef.current.getBoundingClientRect();
    const percentage = Math.max(0, Math.min(100, ((e.clientX - rect.left) / rect.width) * 100));
    const newValue = getValueFromPercentage(percentage);

    if (isDragging === 'min') {
      const newMin = Math.min(newValue, maxValue);
      if (onRangeChange) {
        onRangeChange(newMin, maxValue);
      }
    } else if (isDragging === 'max') {
      const newMax = Math.max(newValue, minValue);
      if (onRangeChange) {
        onRangeChange(minValue, newMax);
      }
    }
  };

  const handleMouseUp = () => {
    setIsDragging(null);
  };

  useEffect(() => {
    if (isDragging) {
      document.addEventListener('mousemove', handleMouseMove);
      document.addEventListener('mouseup', handleMouseUp);
      return () => {
        document.removeEventListener('mousemove', handleMouseMove);
        document.removeEventListener('mouseup', handleMouseUp);
      };
    }
  }, [isDragging, minValue, maxValue, actualMin, actualMax]);

  // Format the numbers for display
  const formatNumber = (num) => {
    if (num === null || num === undefined) return '';
    
    // Use scientific notation for very large or very small numbers
    if (Math.abs(num) >= 1000000 || (Math.abs(num) < 0.001 && num !== 0)) {
      return num.toExponential(2);
    }
    
    // Use fixed notation with appropriate precision
    const precision = Math.abs(num) >= 1 ? 2 : 4;
    return num.toFixed(precision);
  };

  const minPercentage = minValue !== null ? getPercentage(minValue) : 0;
  const maxPercentage = maxValue !== null ? getPercentage(maxValue) : 100;
  const isFiltered = minValue > actualMin || maxValue < actualMax;

  const getMinLabelStyle = () => {
    return {
      left: `${minPercentage}%`,
      top: '20px',
      transform: 'translateX(-100%)' // Always align to the left/outside
    };
  };
  
  const getMaxLabelStyle = () => {
    return {
      left: `${maxPercentage}%`,
      top: '20px',
      transform: 'translateX(0)' // Always align to the right/outside
    };
  };

  return (
    <div className="flex flex-col items-left my-2 justify-between">
      <div className="flex flex-row items-center justify-between mb-2">
        <label className="text-sm text-gray-600 font-medium">
          Filter points by {title}
        </label>
        <button
          onClick={() => onRangeChange && onRangeChange(actualMin, actualMax)}
          className="text-xs text-gray-400 hover:text-gray-600 underline"
          title="Reset filters"
        >
          Reset
        </button>
      </div>
      
      <div className="flex flex-col space-y-3">
        {/* Range Slider */}
        <div className="relative pb-6 px-2">
          <div 
            ref={sliderRef}
            className="relative h-2 bg-gray-200 rounded-lg cursor-pointer mt-3 mb-2 mx-2"
          >
            {/* Active range track */}
            <div
              className={`absolute h-2 rounded-lg ${isFiltered ? 'bg-blue-600' : 'bg-gray-200'}`}
              style={{
                left: `${minPercentage}%`,
                width: `${maxPercentage - minPercentage}%`
              }}
            />
            
            {/* Min handle */}
            <div
              className="absolute w-4 h-4 bg-blue-600 rounded-full cursor-grab active:cursor-grabbing transform -translate-y-1 -translate-x-2 transition-colors shadow-md hover:bg-blue-700"
              style={{ left: `${minPercentage}%` }}
              onMouseDown={handleMouseDown('min')}
              title={`Min: ${formatNumber(minValue)}`}
            />
            
            {/* Min value label */}
            <div
              className="absolute text-xs text-gray-500"
              style={getMinLabelStyle()}
            >
              {formatNumber(minValue)}
            </div>
            
            {/* Max handle */}
            <div
              className="absolute w-4 h-4 bg-blue-600 rounded-full cursor-grab active:cursor-grabbing transform -translate-y-1 -translate-x-2 transition-colors shadow-md hover:bg-blue-700"
              style={{ left: `${maxPercentage}%` }}
              onMouseDown={handleMouseDown('max')}
              title={`Max: ${formatNumber(maxValue)}`}
            />
            
            {/* Max value label */}
            <div
              className="absolute text-xs text-gray-500"
              style={getMaxLabelStyle()}
            >
              {formatNumber(maxValue)}
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};
