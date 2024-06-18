
export default function Checkbox({ text, id, onChange, checked, children }) {
    return (
        <div className="inline-flex items-center">
            <label className="relative flex items-center px-3 rounded-full cursor-pointer" htmlFor={id}>
                <input
                    type="checkbox"
                    className="before:content[''] peer relative h-5 w-5 cursor-pointer appearance-none 
                    rounded-md border border-slate-300 transition-all 
                    before:absolute before:top-2/4 before:left-2/4 before:block before:h-10 
                    before:w-10 before:-translate-y-2/4 before:-translate-x-2/4 
                    before:rounded-full before:bg-blue-gray-500 before:opacity-0 
                    before:transition-opacity checked:border-blue-600 checked:bg-blue-600 
                    checked:before:bg-blue-600 hover:before:opacity-10"
                    id={id}
                    checked={checked}
                    onChange={(e) => onChange(e.target.checked) } 
                    />
                <span
                    className="absolute text-white transition-opacity opacity-0 pointer-events-none top-2/4 left-2/4 -translate-y-2/4 -translate-x-2/4 peer-checked:opacity-100">
                    <svg xmlns="http://www.w3.org/2000/svg" className="h-3.5 w-3.5" viewBox="0 0 20 20" fill="currentColor"
                        stroke="currentColor" strokeWidth="1">
                        <path fillRule="evenodd"
                            d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
                            clipRule="evenodd"></path>
                    </svg>
                </span>
            </label>
            <label className="mt-px text-sm text-gray-500 cursor-pointer select-none" htmlFor={id}>
                {text}
            </label>
        </div>
    )
}


export function Radio({ text, id, onChange, children }) {
    return (
        <div className="flex gap-10">

            <div className="inline-flex items-center">
                <label className="relative flex items-center p-3 rounded-full cursor-pointer" htmlFor={id}>
                    <input name="type"
                        type="radio"
                        className="before:content[''] peer relative h-5 w-5 cursor-pointer 
                    appearance-none rounded-full border border-blue-gray-200 
                    text-gray-900 transition-all before:absolute before:top-2/4 
                    before:left-2/4 before:block before:h-10 before:w-10 
                    before:-translate-y-2/4 before:-translate-x-2/4 
                    before:rounded-full before:bg-blue-gray-500 before:opacity-0 
                    before:transition-opacity checked:border-gray-900 
                    checked:before:bg-gray-900 hover:before:opacity-10"
                        id={id}
                        onChange={onChange()} />
                    <span
                        className="absolute text-gray-900 transition-opacity opacity-0 pointer-events-none 
                    top-2/4 left-2/4 -translate-y-2/4 -translate-x-2/4 peer-checked:opacity-100">
                        <svg xmlns="http://www.w3.org/2000/svg" className="h-3.5 w-3.5" viewBox="0 0 16 16" fill="currentColor">
                            <circle data-name="ellipse" cx="8" cy="8" r="8"></circle>
                        </svg>
                    </span>
                </label>
                <label className="mt-px font-light text-gray-700 cursor-pointer select-none" htmlFor={id}>
                    {text}
                </label>
            </div>
            <div className="inline-flex items-center">
                <label className="relative flex items-center p-3 rounded-full cursor-pointer" htmlFor={id}>
                    <input name="type"
                        type="radio"
                        className="before:content[''] peer relative h-5 w-5 cursor-pointer 
                appearance-none rounded-full border border-blue-gray-200 
                text-gray-900 transition-all before:absolute before:top-2/4 
                before:left-2/4 before:block before:h-10 before:w-10 
                before:-translate-y-2/4 before:-translate-x-2/4 
                before:rounded-full before:bg-blue-gray-500 before:opacity-0 
                before:transition-opacity checked:border-gray-900 
                checked:before:bg-gray-900 hover:before:opacity-10"
                        id={id}
                        onChange={onChange()} />
                    <span
                        className="absolute text-gray-900 transition-opacity opacity-0 pointer-events-none 
                top-2/4 left-2/4 -translate-y-2/4 -translate-x-2/4 peer-checked:opacity-100">
                        <svg xmlns="http://www.w3.org/2000/svg" className="h-3.5 w-3.5" viewBox="0 0 16 16" fill="currentColor">
                            <circle data-name="ellipse" cx="8" cy="8" r="8"></circle>
                        </svg>
                    </span>
                </label>
                <label className="mt-px font-light text-gray-700 cursor-pointer select-none" htmlFor={id}>
                    {text}
                </label>
            </div>
        </div>
    )
}

