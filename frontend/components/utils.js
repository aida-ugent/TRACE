import { Component, useState } from 'react';
import Select from 'react-select'

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
export function saveAsPng(scatterplot, filename='scatter.png') {
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