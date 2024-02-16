import { useState, useLayoutEffect, useRef } from "react";

export function DefaultButton(props) {

  const { onClick, children } = props;

  return (
    <button
      type="button"
      className="min-w-[100px] h-fit select-none justify-center items-center 
      rounded-md bg-slate-500 text-gray-50 mx-1 my-1 px-4 pb-1 pt-1.5 text-sm font-medium leading-normal
       border-0 drop-shadow-sm border border-gray-600 hover:no-underline hover:bg-slate-800
       hover:opacity-75 hover:text-gray-50 focus:opacity-100 focus:shadow-none 
      focus:outline-none disabled:opacity-75 disabled:bg-gray-700 disabled:text-gray-50"
      onClick={onClick}
    >
      {children}
    </button>
  )
}


export function AsyncButton(props) {

  const { onClick, children } = props;
  const [isLoading, setLoading] = useState(false);
  const [buttonWidth, setWidth] = useState("0");
  const ref = useRef();

  const test_dimensions = () => {
    if (ref.current) {
      setWidth((ref.current.offsetWidth + 1).toString());
    }
  }

  useLayoutEffect(() => {
    test_dimensions();
  }, []);

  const handleClick = () => {
    setLoading(true);
    onClick()
      .then(res => setLoading(false));
  }

  return (
    <>
      <button
        type="button"
        className="min-w-[100px] h-fit select-none justify-center items-center 
        rounded-md bg-slate-500 text-gray-50 mx-1 my-1 px-4 pb-1 pt-1.5 text-sm font-medium leading-normal
         border-0 drop-shadow-sm border border-gray-600 hover:no-underline hover:bg-slate-800
         hover:opacity-75 hover:text-gray-50 focus:opacity-100 focus:shadow-none 
        focus:outline-none disabled:opacity-75 disabled:bg-gray-700 disabled:text-gray-50"
        onClick={handleClick}
        disabled={isLoading}
        ref={ref}
      >
        {isLoading ? "loading..." : children}
      </button>
    </>
  )

}

export function ChevronButton(props) {
  const { innerText, onChange, ...other } = props;
  const [visibility, setVisibility] = useState("hidden")

  const toggleVisibility = () => {
    if (visibility == "visible") {
      setVisibility("hidden")
    }
    else {
      setVisibility("visible")
    }
  }

  return (
    <>
      <DefaultButton onClick={toggleVisibility}>
        {<>
          <p className="select-none">{innerText}</p>
          {visibility == "visible" ? (<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-4 h-4 ml-2">
            <path fillRule="evenodd" d="M11.47 7.72a.75.75 0 011.06 0l7.5 7.5a.75.75 0 11-1.06 1.06L12 9.31l-6.97 6.97a.75.75 0 01-1.06-1.06l7.5-7.5z" clipRule="evenodd" />
          </svg>) : (<svg xmlns="http://www.w3.org/2000/svg"
            viewBox="0 0 24 24" fill="currentColor" className="w-4 h-4 ml-2">
            <path fillRule="evenodd" d="M12.53 16.28a.75.75 0 01-1.06 0l-7.5-7.5a.75.75 0 011.06-1.06L12 14.69l6.97-6.97a.75.75 0 111.06 1.06l-7.5 7.5z" clipRule="evenodd" />
          </svg>)}
        </>
        }
      </DefaultButton>
    </>
  )
}

export function ResetButton(props) {
  return (
    <button
      type="button"
      className="select-none inline-block rounded-full mx-2 bg-white/80 p-2 text-sm font-medium leading-normal 
      text-white border-none hover:no-underline focus:opacity-100 focus:shadow-none focus:outline-none"
      onClick={props.onClick}
    >
      {/* <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-6 h-6">
        <path fillRule="evenodd" d="M3.22 3.22a.75.75 0 011.06 0l3.97 3.97V4.5a.75.75 0 011.5 0V9a.75.75 0 01-.75.75H4.5a.75.75 0 010-1.5h2.69L3.22 4.28a.75.75 0 010-1.06zm17.56 0a.75.75 0 010 1.06l-3.97 3.97h2.69a.75.75 0 010 1.5H15a.75.75 0 01-.75-.75V4.5a.75.75 0 011.5 0v2.69l3.97-3.97a.75.75 0 011.06 0zM3.75 15a.75.75 0 01.75-.75H9a.75.75 0 01.75.75v4.5a.75.75 0 01-1.5 0v-2.69l-3.97 3.97a.75.75 0 01-1.06-1.06l3.97-3.97H4.5a.75.75 0 01-.75-.75zm10.5 0a.75.75 0 01.75-.75h4.5a.75.75 0 010 1.5h-2.69l3.97 3.97a.75.75 0 11-1.06 1.06l-3.97-3.97v2.69a.75.75 0 01-1.5 0V15z" clipRule="evenodd" />
      </svg> */}
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-8 h-8 fill-gray-500 hover:fill-gray-800">
        <path d="M6 3a3 3 0 0 0-3 3v1.5a.75.75 0 0 0 1.5 0V6A1.5 1.5 0 0 1 6 4.5h1.5a.75.75 0 0 0 0-1.5H6ZM16.5 3a.75.75 0 0 0 0 1.5H18A1.5 1.5 0 0 1 19.5 6v1.5a.75.75 0 0 0 1.5 0V6a3 3 0 0 0-3-3h-1.5ZM12 8.25a3.75 3.75 0 1 0 0 7.5 3.75 3.75 0 0 0 0-7.5ZM4.5 16.5a.75.75 0 0 0-1.5 0V18a3 3 0 0 0 3 3h1.5a.75.75 0 0 0 0-1.5H6A1.5 1.5 0 0 1 4.5 18v-1.5ZM21 16.5a.75.75 0 0 0-1.5 0V18a1.5 1.5 0 0 1-1.5 1.5h-1.5a.75.75 0 0 0 0 1.5H18a3 3 0 0 0 3-3v-1.5Z" />
      </svg>
    </button>
  )
}

export function ZoomButton(props) {
  return (
    <button
      type="button"
      className="select-none inline-block rounded-full bg-gray-600 mx-2 p-2 text-sm font-medium leading-normal text-white border-none hover:no-underline hover:opacity-75 focus:opacity-100 focus:shadow-none focus:outline-none"
      onClick={props.onClick}
    >
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-6 h-6">
        <path fillRule="evenodd" d="M10.5 3.75a6.75 6.75 0 100 13.5 6.75 6.75 0 000-13.5zM2.25 10.5a8.25 8.25 0 1114.59 5.28l4.69 4.69a.75.75 0 11-1.06 1.06l-4.69-4.69A8.25 8.25 0 012.25 10.5z" clipRule="evenodd" />
      </svg>
    </button>
  )
}

export function SettingsButton(props) {
  return (
    <button
      type="button"
      className="select-none inline-block rounded-full bg-slate-500 mx-2 p-2 text-sm font-medium 
      leading-normal text-white border-none hover:no-underline hover:opacity-75 focus:opacity-100 
      focus:shadow-none focus:outline-none hover:bg-slate-800"
      onClick={props.onClick}
    >
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-6 h-6">
        <path d="M18.75 12.75h1.5a.75.75 0 000-1.5h-1.5a.75.75 0 000 1.5zM12 6a.75.75 0 01.75-.75h7.5a.75.75 0 010 1.5h-7.5A.75.75 0 0112 6zM12 18a.75.75 0 01.75-.75h7.5a.75.75 0 010 1.5h-7.5A.75.75 0 0112 18zM3.75 6.75h1.5a.75.75 0 100-1.5h-1.5a.75.75 0 000 1.5zM5.25 18.75h-1.5a.75.75 0 010-1.5h1.5a.75.75 0 010 1.5zM3 12a.75.75 0 01.75-.75h7.5a.75.75 0 010 1.5h-7.5A.75.75 0 013 12zM9 3.75a2.25 2.25 0 100 4.5 2.25 2.25 0 000-4.5zM12.75 12a2.25 2.25 0 114.5 0 2.25 2.25 0 01-4.5 0zM9 15.75a2.25 2.25 0 100 4.5 2.25 2.25 0 000-4.5z" />
      </svg>
    </button>
  )
}


export function ChevronRightButton(props) {
  return (
    <button
      type="button"
      className="select-none inline-block rounded-full bg-slate-500 p-2 text-sm font-medium 
      leading-normal text-white border-none hover:no-underline hover:opacity-75 focus:opacity-100 
      focus:shadow-none focus:outline-none hover:bg-slate-800"
      onClick={props.onClick}
    >
      <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" className="w-6 h-6">
        <path fillRule="evenodd" d="M4.72 3.97a.75.75 0 011.06 0l7.5 7.5a.75.75 0 010 1.06l-7.5 7.5a.75.75 0 01-1.06-1.06L11.69 12 4.72 5.03a.75.75 0 010-1.06zm6 0a.75.75 0 011.06 0l7.5 7.5a.75.75 0 010 1.06l-7.5 7.5a.75.75 0 11-1.06-1.06L17.69 12l-6.97-6.97a.75.75 0 010-1.06z" clipRule="evenodd" />
      </svg>
    </button>
  )
}