import { Component, useState } from 'react';
import Select from 'react-select'
import { SubmitButton } from './buttons';

export function getPointSize(numPoints) {
  // point size depending on number of points
  // between [0.1, 10]
  return Math.min(10, Math.max(0.1, 20 * Math.exp(-0.458145 * Math.log10(numPoints))));
}


export function prettyPrint(value) {
  if (typeof value === 'string') {
    return value;
  } else if (typeof value === 'number') {
    return value.toPrecision(3);
  } else {
    return value;
  }
}

// from https://github.com/flekschas/regl-scatterplot/blob/39d353d5cf0f0e37f821c498322773989a1f5d1d/example/utils.js#L19
export function downloadBlob(blob, name = 'file.txt') {
  const link = document.createElement('a');
  link.href = URL.createObjectURL(blob);
  link.download = name;

  document.body.appendChild(link);

  link.dispatchEvent(
    new MouseEvent('click', {
      bubbles: true,
      cancelable: true,
      view: window,
    })
  );

  document.body.removeChild(link);
}

function getCanvasWithBackground(scatterplot, bgColor) {
  const canvas = scatterplot.get('canvas');
  var width = canvas.width;
  var height = canvas.height;

  const combinedCanvas = document.createElement("canvas");
  combinedCanvas.width = width;
  combinedCanvas.height = height;

  const combinedCtx = combinedCanvas.getContext('2d');
  combinedCtx.fillStyle = bgColor;
  combinedCtx.fillRect(0, 0, width, height);
  combinedCtx.drawImage(canvas, 0, 0, width, height);
  return combinedCanvas;
}

// from https://github.com/flekschas/regl-scatterplot/blob/39d353d5cf0f0e37f821c498322773989a1f5d1d/example/utils.js#L19
export function saveAsPng(scatterplot, filename = 'scatter.png') {
  const imageObject = new Image();
  imageObject.onload = () => {
    getCanvasWithBackground(scatterplot, '#ffffff').toBlob((blob) => {
      downloadBlob(blob, filename);
    });
  };
  imageObject.src = scatterplot.get('canvas').toDataURL('image/png', 1.0);
}



class NamedSlider extends Component {
  constructor(props) {
    super(props);
    this.state = {
      value: props.defaultValue,
    };

    this.showCurrentValue = (value) => {
      this.setState({ value: value });
      this.props.onChange(value)
    };
  }

  render() {
    return (
      <>
        <label
          className="text-sm text-gray-500 w-fit min-w-fit"
          htmlFor="neighborSlider">
          {
            this.props.names !== undefined ?
              this.props.label + " " + this.props.names[this.state.value] :
              this.props.label + " " + this.state.value
          }
        </label>
        <input
          id="neighborSlider"
          className="transparent h-[2px] cursor-pointer appearance-none border-transparent bg-neutral-300 my-2"
          type="range" min={this.props.min} max={this.props.max} step={this.props.step} defaultValue={this.props.defaultValue}
          onChange={(event) => this.showCurrentValue(+event.target.value)} />
      </>
    );
  }
}

export default NamedSlider

export const Slider = ({ onChange, min, max, step, defaultValue, disabled, id }) => {
  return (
    <input
      className={`transparent h-[2px] appearance-none border-transparent bg-neutral-300 
    mb-2 mt-3 ${disabled ? 'accent-slate-100' : 'cursor-pointer'}`}
      type="range"
      min={min}
      max={max}
      step={step}
      defaultValue={defaultValue}
      onChange={(event) => onChange(+event.target.value)}
      id={id}
      disabled={disabled}
    />
  )
}

const RadioButton = ({ onChange, value, checked }) => (
  <label>
    <input type="radio" name="radio-button-group" value={value} onChange={onChange} checked={checked} /> {value}
  </label>
);


export function RadioSelect(props) {
  const { onChange, defaultValue, options, ...other } = props;
  const [currentValue, setCurrentValue] = useState(defaultValue);

  function onRadioChange(event) {
    setCurrentValue(event.target.value);
    onChange(event.target.value);
  }

  return (
    <div className='flex flex-col w-fit h-fit'>
      {
        options.map(option =>
          <RadioButton key={option} value={option} checked={option == currentValue} onChange={onRadioChange} />
        )
      }
    </div>
  )
}


export function DropdownSelectOld(props) {
  const { options, selected, onChange, id, ...other } = props

  return (
    <select id={id} value={selected} onChange={(event) => onChange(event.target.value)}>
      {options.map((option, optionIdx) => (
        <option key={optionIdx} value={option}>
          {option}
        </option>
      ))}
    </select>
  );
}


export function ReactSelect({ options, selected, onChange, isDisabled = false, menuPlacement = 'auto' }) {

  return (
    <Select
      className="min-w-max text-slate-600 text-left"
      options={options}
      isClearable={false}
      isSearchable={true}
      isDisabled={isDisabled}
      menuPlacement={menuPlacement}
      value={{ 'value': selected, 'label': selected }}
      onChange={(selection) => { onChange(selection["value"]) }}
    />
  )
}



export function HoverNote(props) {
  const { visible, color, position, children } = props;
  return (
    <button
      type="button"
      className={"select-none fixed rounded-md bg-white/90 outline outline-2 text-left m-2 p-2 text-base font-medium leading-normal text-black"}
      style={{
        'outlineColor': color,
        'visibility': visible,
        'top': position[1],
        'left': position[0],
      }}
    >
      {children}
    </button >
  )
}



export function SavePointForm(props) {
  const { scatterplot, children } = props;
  const [name, setName] = useState('selection name');

  const handleSubmit = (e) => {
    e.preventDefault();
    console.log(`Form submitted, ${name}`);

    const selectedPoints = scatterplot.get('selectedPoints');

    if (selectedPoints.length > 1) {
      fetch("/backend/savePointSelection", {
        method: "POST",
        body: JSON.stringify({
          points: selectedPoints,
          selection_name: name
        }),
        headers: {
          "Content-type": "application/json; charset=UTF-8"
        }
      })
        .then(response => response.json())
        .then(data => {
          console.log('Success:', data);
        })
        .catch((error) => {
          console.error('Error:', error);
        });
    }
  }

  return (
    <div className='flex flex-wrap items-start justify-between mb-2'>
      <form onSubmit={handleSubmit} className='flex flex-wrap justify-between items-start'>
        <input
          onChange={(e) => setName(e.target.value)}
          value={name}
          type="text"
          name="name"
          className="rounded-md bg-white  px-2 pb-2 pt-2.5 w-3/5 text-sm font-medium leading-normal text-gray-700 "
          minLength={1} />
        <SubmitButton>
          save
        </SubmitButton>
      </form>
    </div>
  );

}